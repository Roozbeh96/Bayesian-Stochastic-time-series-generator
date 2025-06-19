x_test = torch.rand(1, 1)*0-5.0
pi, mu, sigma = Given_min_model(x_test)


# Sample 'num_samples' values for each of the 500 x_test values
num_samples = 50  # Input parameter for number of samples to generate
y_samples = []

for i in range(x_test.shape[0]):
    y_samples.append(sample_from_combined_pdf(pi[i:i+1, :], mu[i:i+1, :], sigma[i:i+1, :], num_samples).detach().numpy())

# Plotting the results
plt.figure(figsize=(8, 6))
plt.scatter(x_train.numpy(), y_train.numpy(), alpha=0.5, label='Training Data')

# Loop over x_test and corresponding y_samples to plot individual samples
for i, y_s in enumerate(y_samples):
    plt.scatter(np.repeat(x_test[i].numpy(), num_samples), y_s, s=2, label=f'Sample {i+1}')

plt.title('Mixture Density Network: Sampled Outputs')
plt.xlabel('x')
plt.ylabel('y')
plt.show()