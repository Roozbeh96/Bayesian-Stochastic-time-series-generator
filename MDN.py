import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from torch.distributions import Normal, Categorical
import scipy.io as sio
import h5py
from torch.distributions import Normal, LogNormal

# For a sample problem
 
# def data_generator():
#     t = np.linspace(0, 1, 1000)  # from 0 to 1 seconds, 1000 points
#     # Add noise (Gaussian)
#     noise = np.random.normal(0, 0.05, size=t.shape)  # mean=0, std=0.05

#     # Compute x using the formula
#     x = t + 0.3 * np.sin(2 * np.pi * t) + noise
#     return x.reshape(-1,1), t.reshape(-1,1)

# x_train, y_train = data_generator()

# # Plot the data with shaded std deviation
# plt.figure(figsize=(8, 5))
# plt.scatter(x_train, y_train, alpha=0.5, label="Noisy Data")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.legend()
# plt.title("Inverse Problem")
# plt.show()

# x_train = torch.tensor(x_train, dtype=torch.float32)
# y_train = torch.tensor(y_train, dtype=torch.float32)

# Define the Mixture Density Network
class MDN(nn.Module):
    def __init__(self, input_dim=1, output_dim=1, num_gaussians=5):
        super(MDN, self).__init__()
        self.num_gaussians = num_gaussians
        self.fc = nn.Sequential(
            nn.Linear(input_dim, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU(),
            nn.Linear(32, 32),
            nn.ReLU()
        )
        self.pi_layer = nn.Linear(32, num_gaussians)  # Mixing coefficients
        self.mu_layer = nn.Linear(32, num_gaussians)  # Means of Gaussians
        self.sigma_layer = nn.Linear(32, num_gaussians)  # Standard deviations
    
    def forward(self, x):
        h = self.fc(x)
        pi = torch.softmax(self.pi_layer(h), dim=1)  # Mixing coefficients (softmax)
        mu = self.mu_layer(h)  # Means
        sigma =  torch.exp(self.sigma_layer(h)).clamp(min=0.005)  # Standard deviations (exp for positivity)
        return pi, mu, sigma

# Define the Negative Log Likelihood loss for MDN
def mdn_loss(pi, mu, sigma, y, disttype = "Normal"):
    y = y.expand(-1, mu.shape[-1])  # Fix dimension expansion
    if disttype == "Normal":
        normal_dist = Normal(mu, sigma)
        likelihoods = torch.exp(normal_dist.log_prob(y))  # Due to numerical instability we first take a log and then exp
    else:
        log_normal_dist = LogNormal(mu, sigma)
        likelihoods = torch.exp(log_normal_dist.log_prob(y))
    weighted_likelihoods = pi * likelihoods  # Weight by mixing coefficients
    loss = -torch.log(weighted_likelihoods.sum(dim=1) + 1e-8).mean()  # Log sum exp trick
    return loss

# Training the model
# sample_model = MDN()
# optimizer = optim.Adam(sample_model.parameters(), lr=0.01)

# # Training loop
# epochs = 10000
# for epoch in range(epochs):
#     optimizer.zero_grad()
#     pi, mu, sigma = sample_model(x_train)
#     dist_type = "Normal"
#     loss = mdn_loss(pi, mu, sigma, y_train, dist_type)
#     loss.backward()
#     optimizer.step()
#     if epoch % 50 == 0:
#         print(f'Epoch {epoch}, Loss: {loss.item():.4f}')

# Final PDF calculation function (same as before)
def final_pdf(pi, mu, sigma, x, disttype="Normal"):
    pdf_values = torch.zeros(x.shape[0], pi.shape[1])  # [1000, 5]
    
    # Compute the individual PDF values for each Gaussian component
    for i in range(pi.shape[1]):  # num_gaussians
        if disttype=="Normal":
            dist = Normal(mu[:, i], sigma[:, i])
            pdf_values[:, i] = dist.log_prob(x).exp()  # Exponentiate to get the actual PDF
        else:
            dist = LogNormal(mu[:, i], sigma[:, i])
            pdf_values[:, i] = dist.log_prob(x).exp()  # Exponentiate to get the actual PDF
    
    # Weight the PDFs by the mixing coefficients (pi)
    weighted_pdf = pdf_values * pi  # Broadcasting pi across [1000, 5]
    
    # Sum the weighted PDFs over all components
    final_pdf_values = weighted_pdf.sum(dim=1)  
    # for Nomral (final_pdf_values*(x[1]-x[0])).sum() must be 1
    # for LogNormal torch.trapz(final_pdf_values, x)
    return final_pdf_values

# Sampling from the combined PDF based on the learned model parameters
def sample_from_combined_pdf(pi, mu, sigma, num_samples=100, disttype="Normal"):
    # print(num_samples)
    """
    Sample 'num_samples' values from the final combined distribution using the CDF.
    
    pi: mixing coefficients (5 values)
    mu: means of the 5 Gaussians
    sigma: standard deviations of the 5 Gaussians
    num_samples: number of samples to generate (default is 100)
    """
    if disttype=="Normal":
        x_values = torch.linspace(-10, 10, steps=1000)  # Define the range for sampling
    else:
        x_values = torch.logspace(-5, 5, steps=1000)  # Define the range for sampling

    
    pdf_values = final_pdf(pi, mu, sigma, x_values, disttype)  # Compute the PDF
    
    # Normalize the PDF to get the CDF
    cdf_values = torch.cumsum(pdf_values, dim=0)
    cdf_values = cdf_values / cdf_values[-1]  # Normalize so that CDF ends at 1
    
    # Generate 'num_samples' random values to sample from the CDF
    u = torch.rand(num_samples)  # Generate random values to sample from the CDF
    
    samples = []
    for val in u:
        sample_index = (cdf_values > val).nonzero(as_tuple=True)[0][0].item()
        sample = x_values[sample_index]
        samples.append(sample)
    
    return torch.tensor(samples)

# Example usage:
# x_test = torch.rand(500, 1)  # 500 values of x_test (batch of 500)
# pi, mu, sigma = sample_model(x_test)  # Get the corresponding pi, mu, sigma

# # Sample 'num_samples' values for each of the 500 x_test values
# num_samples = 50  # Input parameter for number of samples to generate
# y_samples = []

# for i in range(x_test.shape[0]):
#     dist_type = "Normal"
#     y_samples.append(sample_from_combined_pdf(pi[i:i+1, :], mu[i:i+1, :], sigma[i:i+1, :], num_samples,disttype = dist_type).detach().numpy())

# # Plotting the results
# plt.figure(figsize=(8, 6))
# plt.scatter(x_train.numpy(), y_train.numpy(), alpha=0.5, label='Training Data')

# # Loop over x_test and corresponding y_samples to plot individual samples
# for i, y_s in enumerate(y_samples):
#     plt.scatter(np.repeat(x_test[i].numpy(), num_samples), y_s, s=2, label=f'Sample {i+1}')

# plt.title('Mixture Density Network: Sampled Outputs')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()


# Given min value, distribution of max

HotWireWT7 = h5py.File(r'SGFV_param_HotWT7.mat','r')
T_domain_WT7 = np.array(HotWireWT7.get('T_domain_WT7')).T
X_domain_WT7 = np.array(HotWireWT7.get('X_domain_WT7')).T
res_WT7 = np.array(HotWireWT7.get('res_WT7'))

VF_HotWT7 = HotWireWT7['VF_HotWT7']
VF_HotWT7_u_tau = VF_HotWT7['u_tau'][0]


duNN_refs = HotWireWT7['duNN_HotWT7']
maxmin_pairs_cell_HotWT7_refs = HotWireWT7['value_pairs_cell_HotWT7']
duwrtlocalmaxmin_HotWT7_refs = HotWireWT7['duwrtlocalmaxmin_HotWT7']

# To dereference:
duNN_HotWT7_list = []

for ref in duNN_refs:
    # Iterate over each element inside 'ref'
    for i, subref in enumerate(ref):
        # subref is now a real HDF5 reference
        dereferenced_data = HotWireWT7[subref]
        array_data = np.array(dereferenced_data)
        duNN_HotWT7_list.append(array_data)

duwrtlocalmaxmin_HotWT7_list = []

for ref in duwrtlocalmaxmin_HotWT7_refs:
    # Iterate over each element inside 'ref'
    for i, subref in enumerate(ref):
        # subref is now a real HDF5 reference
        dereferenced_data = HotWireWT7[subref]
        array_data = np.array(dereferenced_data)
        duwrtlocalmaxmin_HotWT7_list.append(array_data)

maxmin_pairs_HotWT7_list = []

for ref in maxmin_pairs_cell_HotWT7_refs:
    # Iterate over each element inside 'ref'
    for i, subref in enumerate(ref):
        # subref is now a real HDF5 reference
        dereferenced_data = HotWireWT7[subref]
        array_data = np.array(dereferenced_data)
        maxmin_pairs_HotWT7_list.append(array_data)


min_values = maxmin_pairs_HotWT7_list[0][1]/VF_HotWT7_u_tau
max_values = maxmin_pairs_HotWT7_list[0][0]/VF_HotWT7_u_tau


# #Training the model
# Given_min_model = MDN()
# optimizer = optim.Adam(Given_min_model.parameters(), lr=0.01)

# x_train = torch.tensor(min_values.reshape(-1,1), dtype=torch.float32)
# y_train = torch.tensor(max_values.reshape(-1,1), dtype=torch.float32)

# # Training loop
# epochs = 10000
# for epoch in range(epochs):
#     optimizer.zero_grad()
#     pi, mu, sigma = Given_min_model(x_train)
#     dist_type = "Normal"
#     loss = mdn_loss(pi, mu, sigma, y_train, dist_type)
#     loss.backward()
#     optimizer.step()
#     if epoch % 50 == 0:
#         print(f'Epoch {epoch}, Loss: {loss.item():.4f}')

'''If we want do inferencing'''
# Given_min_model.load_state_dict(torch.load('G:\My Drive\Research\VFfeaturedVorX\Given_min_model.pth'))

# Given_min_model.eval()

# x_test = torch.rand(1, 1)*0-4.0
# pi, mu, sigma = Given_min_model(x_test)


# # Sample 'num_samples' values for each of the 500 x_test values
# num_samples = 50  # Input parameter for number of samples to generate
# y_samples = []

# for i in range(x_test.shape[0]):
#     dist_type = "Normal"
#     y_samples.append(sample_from_combined_pdf(pi[i:i+1, :], mu[i:i+1, :], sigma[i:i+1, :], num_samples,disttype = dist_type).detach().numpy())

# # Plotting the results
# plt.figure(figsize=(8, 6))
# plt.scatter(x_train.numpy(), y_train.numpy(), alpha=0.5, label='Training Data')

# # Loop over x_test and corresponding y_samples to plot individual samples
# for i, y_s in enumerate(y_samples):
#     plt.scatter(np.repeat(x_test[i].numpy(), num_samples), y_s, s=2, label=f'Sample {i+1}')

# plt.title('Mixture Density Network: Sampled Outputs')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# # Saving the trained model
# torch.save(Given_min_model.state_dict(), 'Given_min_model.pth')


# Given_max_model = MDN()
# optimizer = optim.Adam(Given_max_model.parameters(), lr=0.01)

# x_train = torch.tensor(max_values.reshape(-1,1), dtype=torch.float32)
# y_train = torch.tensor(min_values.reshape(-1,1), dtype=torch.float32)


# # Training loop
# epochs = 10000
# for epoch in range(epochs):
#     optimizer.zero_grad()
#     pi, mu, sigma = Given_max_model(x_train)
#     dist_type = "Normal"
#     loss = mdn_loss(pi, mu, sigma, y_train, dist_type)
#     loss.backward()
#     optimizer.step()
#     if epoch % 50 == 0:
#         print(f'Epoch {epoch}, Loss: {loss.item():.4f}')

# # Visualizing the learned distribution
'''If we want do inferencing'''
# Given_max_model.load_state_dict(torch.load('G:\My Drive\Research\VFfeaturedVorX\Given_max_model.pth'))

# Given_max_model.eval()
# x_test = torch.rand(1, 1)*0-4.0
# pi, mu, sigma = Given_max_model(x_test)


# # Sample 'num_samples' values for each of the 500 x_test values
# num_samples = 50  # Input parameter for number of samples to generate
# y_samples = []

# for i in range(x_test.shape[0]):
#     dist_type = "Normal"
#     y_samples.append(sample_from_combined_pdf(pi[i:i+1, :], mu[i:i+1, :], sigma[i:i+1, :], num_samples,disttype = dist_type).detach().numpy())

# # Plotting the results
# plt.figure(figsize=(8, 6))
# plt.scatter(y_train.numpy(),x_train.numpy() , alpha=0.5, label='Training Data')

# # Loop over x_test and corresponding y_samples to plot individual samples
# for i, y_s in enumerate(y_samples):
#     plt.scatter(y_s, np.repeat(x_test[i].numpy(), num_samples), s=2, label=f'Sample {i+1}')

# plt.title('Mixture Density Network: Sampled Outputs')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# # Saving the trained model
# torch.save(Given_max_model.state_dict(), 'Given_max_model.pth')

duwrtlocalmaxmin_HotWT7 = duwrtlocalmaxmin_HotWT7_list[0]/VF_HotWT7_u_tau
duNN_HotWT7 = duNN_HotWT7_list[0]/(VF_HotWT7_u_tau*res_WT7[0,0])

Given_dist_NLmin_NLmax_model = MDN()
optimizer = optim.Adam(Given_dist_NLmin_NLmax_model.parameters(), lr=0.01)

x_train = torch.tensor(duwrtlocalmaxmin_HotWT7.reshape(-1,1), dtype=torch.float32)
y_train = torch.log(torch.tensor(duNN_HotWT7.reshape(-1,1), dtype=torch.float32)+1e-5)
# y_train = torch.tensor(duNN_HotWT7.reshape(-1,1), dtype=torch.float32)

# nonzero_mask = (y_train != 0).squeeze()

# # Filter both x_train and y_train using the mask
# x_train = x_train[nonzero_mask]
# y_train= y_train[nonzero_mask]


# Training loop
epochs = 10000
for epoch in range(epochs):
    optimizer.zero_grad()
    pi, mu, sigma = Given_dist_NLmin_NLmax_model(x_train)
    dist_type = "Normal"
    loss = mdn_loss(pi, mu, sigma, y_train, dist_type)
    loss.backward()
    optimizer.step()
    if epoch % 50 == 0:
        print(f'Epoch {epoch}, Loss: {loss.item():.4f}')

'''If we want do inferencing'''
Given_dist_NLmin_NLmax_model.load_state_dict(torch.load('G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.pth'))
Given_dist_NLmin_NLmax_model.eval()

x_test = torch.rand(1, 1)*0+1.0
pi, mu, sigma = Given_dist_NLmin_NLmax_model(x_test)


# Sample 'num_samples' values for each of the 500 x_test values
num_samples = 50  # Input parameter for number of samples to generate
y_samples = []

for i in range(x_test.shape[0]):
    dist_type = "Normal"
    y_samples.append(sample_from_combined_pdf(pi[i:i+1, :], mu[i:i+1, :], sigma[i:i+1, :], num_samples,disttype = dist_type).detach().numpy())

# Plotting the results
plt.figure(figsize=(8, 6))
plt.scatter(x_train.numpy(), np.exp(y_train.numpy()), alpha=0.5, label='Training Data')

# Loop over x_test and corresponding y_samples to plot individual samples
for i, y_s in enumerate(y_samples):
    plt.scatter(np.repeat(x_test[i].numpy(), num_samples), np.exp(y_s), s=2, label=f'Sample {i+1}')

plt.yscale('log')
plt.title('Mixture Density Network: Sampled Outputs')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# Saving the trained model
torch.save(Given_dist_NLmin_NLmax_model.state_dict(), 'Given_dist_NLmin_NLmax_model.pth')