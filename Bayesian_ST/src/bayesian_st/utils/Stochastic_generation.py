import numpy as np
from utils.velocity_flu_gen import velocity_fluctuating_generation

class stochastic_generation:
    def __init__(self, *, u_tau, z_0, ks, delta, lambda_T, 
                 nu, z_min, z_max, res_z, res_x):
        self.u_tau = u_tau
        self.z_0 = z_0  
        self.ks = ks
        self.delta = delta
        self.lambda_T = lambda_T
        self.nu = nu
        self.z = np.arange(z_min, z_max+res_z, res_z)
        self.x = np.arange(0, delta+res_x, res_x)
        
    def mean_velocity_field(self, *, numb_frame):
        kappa = 0.39

        U_profile = self.u_tau / kappa * np.log(self.z / self.z_0)
        self.u = np.tile(U_profile[:, np.newaxis, np.newaxis],
                         (1, np.shape(self.x)[0], numb_frame))
        self.w = np.zeros_like(self.u)
        
    def Bayesian_sto(self):

        velocity_fluctuating_generation(self)   

        pass
    