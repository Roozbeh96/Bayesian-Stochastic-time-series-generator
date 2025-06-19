import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from torch.distributions import Normal, Categorical
import scipy.io as sio
import h5py



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
    
Given_max_model = MDN(1,1,5)
Given_max_model.load_state_dict(torch.load('G:\My Drive\Research\VFfeaturedVorX\Given_max_model.pth'))

Given_max_model.eval()

x_test_tensor = torch.tensor([[x_test]], dtype=torch.float32)

pi, mu, sigma = Given_max_model(x_test_tensor)  # Get the corresponding pi, mu, sigma !!Maybe there is an issue with inferencing  

ReturnList = [pi.detach().numpy(), mu.detach().numpy(), sigma.detach().numpy()]