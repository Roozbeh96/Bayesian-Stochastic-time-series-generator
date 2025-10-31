import numpy as np
from utils.Stochastic_generation import stochastic_generation


# Defining Physical properties. 
'''
u_{tau}[m/s](friction velocity),
z_{0}[m](aerodynamic roughness length),
ks[m](equivalent sand-grain roughness height),
delta[m](boundary layer thickness),
lambda[m](Taylor micro-scale),
nu[m^2/s](kinematic viscosity)
'''
z_0 = 6.2e-5
u_tau = 0.39
kappa = 0.39
AFR = 8.5
ks = z_0 / np.exp(-kappa * AFR)
delta = 0.4
lambda_T = 0.01
nu = 1.57e-5

z_min = 50*nu / u_tau
z_max = 0.25*delta
res_z = 0.4*lambda_T/10
Kr = 32
res_x = lambda_T/Kr

# Generating an Object.


Gen_sample = stochastic_generation(u_tau=u_tau, z_0=z_0, ks=ks,
            delta=delta, lambda_T=lambda_T, nu=nu, z_min=z_min,
              z_max=z_max, res_z=res_z, res_x=res_x)

# Generating the mean velocity field
numb_frame = 36
Gen_sample.mean_velocity_field(numb_frame = numb_frame)

# Generating fluctuating velocity signal using Bayesian Stochastic method
Gen_sample.Bayesian_sto()