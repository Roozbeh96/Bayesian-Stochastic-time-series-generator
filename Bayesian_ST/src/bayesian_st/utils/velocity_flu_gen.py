import h5py
import numpy as np
from src.bayesian_st.utils.MDN import MDN_model as MDN
import os
from scipy.stats import norm
import torch


def velocity_fluctuating_generation(obj):

    data1 = h5py.File(r'std_uprime_chopped.mat','r')
    data2 = h5py.File(r'std_wprime_chopped.mat','r')
    std_uprime_O_u_tau_chopped = np.array(data1.get('std_uprime_O_u_tau_chopped')).T
    std_wprime_O_u_tau_chopped = np.array(data2.get('std_wprime_O_u_tau_chopped')).T
    ''' 
    The correlation coeficient between u' and w' can be estimated using
    the chopped_signal_analysis.m line 80-84
    '''
    Number_of_Gaussians = 5
    rho_uw = -0.3 
    path_ = os.getcwd()
    path_ = os.path.join(path_, "Given_max_model.pth")
    Given_min_model = MDN(input_dim= 1,num_gaussians=Number_of_Gaussians)
    Given_min_model.load_state_dict(torch.load(path_))
    Given_min_model.eval()

    path_ = os.getcwd()
    path_ = os.path.join(path_, "Given_min_model.pth")
    Given_max_model = MDN(input_dim=1,num_gaussians=Number_of_Gaussians)
    Given_max_model.load_state_dict(torch.load(path_))
    Given_max_model.eval()

    path_ = os.getcwd()
    path_ = os.path.join(path_, "Given_dist_NLmin_NLmax_model_Long.pth")  
    Given_dist_NLmin_NLmax_model = MDN(input_dim=1,num_gaussians=Number_of_Gaussians)
    Given_dist_NLmin_NLmax_model.load_state_dict(torch.load(path_))
    Given_dist_NLmin_NLmax_model.eval()
    
    if os.path.exists(file_path):
        print("Loading existing array...")
        data = np.load(file_path)
    else:
        print("File not found — computing and saving...")
        # Your code that generates the array
        data = np.random.randn(100, 100)
        np.save(file_path, data)

    size_minmaxvec = 1e6
    minmaxvec = np.zeros((size_minmaxvec,2))
    minind = 0
    maxind = 0
    # take samples from distribution of u' to start the process(Fig3(b))
    start = np.random.normal(loc = 0, scale = 1)*std_uprime_O_u_tau_chopped[0,0]
    sample = np.random.normal(loc = 0, scale = 1)*std_uprime_O_u_tau_chopped[0,0]
    if start > sample:
        minmaxvec[minind,0] = sample
        minind += 1
    else:
        minmaxvec[maxind,1] = sample
        maxind += 1

    is_min_sampled = start > sample

    for x_ in range(0, 2*size_minmaxvec-2):
        if is_min_sampled:
            x_test_tensor = torch.tensor([[minmaxvec[minind-1,0]]], dtype=torch.float32)
            pi, mu, sigma = Given_min_model(x_test_tensor)
            minmaxvec[maxind,1] = minmaxvec[minind-1,0]
            # make sure max > min
            while minmaxvec[maxind,1] <= minmaxvec[minind-1,0]:
                minmaxvec[maxind,1] = 0
                randomvar = np.random.uniform(0,1)
                for j in range(Number_of_Gaussians):
                    minmaxvec[maxind,1] += pi[1,j]*norm.ppf(randomvar,loc = mu[0,j],scale = sigma[0,j])
            maxind += 1
        else:
            x_test_tensor = torch.tensor([[minmaxvec[maxind-1,1]]], dtype=torch.float32)
            pi, mu, sigma = Given_max_model(x_test_tensor)
            minmaxvec[minind,0] = minmaxvec[minind-1,1]
            # make sure max > min
            while minmaxvec[minind,0] >= minmaxvec[maxind-1,1]:
                minmaxvec[minind,0] = 0
                randomvar = np.random.uniform(0,1)
                for j in range(Number_of_Gaussians):
                    minmaxvec[minind,0] += pi[0,j]*norm.ppf(randomvar,loc = mu[0,j],scale = sigma[0,j])
            minind += 1
        is_min_sampled = not is_min_sampled
    
    minind = 0
    maxind = 0

    SGuprime_O_u_tau = np.zeros((obj.z.shape[0], obj.x.shape[0] * obj.u.shape[2]))

    SGuprime = np.zeros((obj.z.shape[0], obj.x.shape[0] * obj.u.shape[2]))
    SGwprime = np.zeros((obj.z.shape[0], obj.x.shape[0] * obj.u.shape[2]))
    res = obj.x[1] - obj.x[0]
    for r in range(obj.z.shape[0]):
        if r == 0:
            start = np.random.normal(loc = 0, scale = 1)*std_uprime_O_u_tau_chopped[0,0]
            SGuprime_O_u_tau[r,0] = start
            if start >= minmaxvec[0,0]:
                is_min_sampled = True
            else:
                is_min_sampled = False
            first_toward = True
            i = 1
            while i <= SGuprime_O_u_tau.shape[1]:
                if first_toward:
                    if is_min_sampled:
                        duNLmin = np.abs(SGuprime_O_u_tau[r,i-1] - minmaxvec[minind,0])
                        x_test_tensor = torch.tensor([[duNLmin]], dtype=torch.float32)
                        pi, mu, sigma = Given_dist_NLmin_NLmax_model(x_test_tensor) 
                        randomvar = np.random.uniform(0,1)
                        duNN = 0
                        for j in range(Number_of_Gaussians):
                            duNN += pi[0,j]*norm.ppf(randomvar,loc = mu[0,j],scale = sigma[0,j])
                        pass
                        if (SGuprime_O_u_tau[r,i-1]-np.exp(duNN)*res)<minmaxvec[minind,0]:
                            SGuprime_O_u_tau[r,i] = minmaxvec[minind,0]
                            minind += 1
                            i += 1
                            is_min_sampled = not is_min_sampled
                            first_toward = False
                        else:
                            SGuprime_O_u_tau[r,i] = SGuprime_O_u_tau[r,i-1]-np.exp(duNN)*res
                            i += 1
                    else:
                        duNLmax = np.abs(SGuprime_O_u_tau[r,i-1] - minmaxvec[maxind,1])
                        x_test_tensor = torch.tensor([[duNLmax]], dtype=torch.float32)
                        pi, mu, sigma = Given_dist_NLmin_NLmax_model(x_test_tensor)
                        randomvar = np.random.uniform(0,1)
                        duNN = 0
                        for j in range(Number_of_Gaussians):
                            duNN += pi[0,j]*norm.ppf(randomvar,loc = mu[0,j],scale = sigma[0,j])
                        if (SGuprime_O_u_tau[r,i-1]+np.exp(duNN)*res)>minmaxvec[maxind,1]:
                            SGuprime_O_u_tau[r,i] = minmaxvec[maxind,1]
                            maxind += 1
                            i += 1
                            is_min_sampled = not is_min_sampled
                            first_toward = False
                        else:
                            SGuprime_O_u_tau[r,i] = SGuprime_O_u_tau[r,i-1]+np.exp(duNN)*res
                            i += 1
                else:
                    if is_min_sampled:
                        duNLmax = np.abs(SGuprime_O_u_tau[r,i-1]-minmaxvec[maxind-1,1])
                        duNLmin = np.abs(SGuprime_O_u_tau[r,i-1]-minmaxvec[minind,0])
                        x_test_tensor = torch.tensor([[min(duNLmin, duNLmax)]], dtype=torch.float32)
                        pi, mu, sigma = Given_dist_NLmin_NLmax_model(x_test_tensor) 
                        randomvar = np.random.uniform(0,1)
                        duNN = 0
                        for j in range(Number_of_Gaussians):
                            duNN += pi[0,j]*norm.ppf(randomvar,loc = mu[0,j],scale = sigma[0,j])
                        if (SGuprime_O_u_tau[r,i-1]-np.exp(duNN)*res)<minmaxvec[minind,0]:
                            SGuprime_O_u_tau[r,i] = minmaxvec[minind,0]
                            minind += 1
                            i += 1
                            is_min_sampled = not is_min_sampled
                        else:
                            SGuprime_O_u_tau[r,i] = SGuprime_O_u_tau[r,i-1]-np.exp(duNN)*res
                            i += 1
                    else:
                        duNLmax = np.abs(SGuprime_O_u_tau[r,i-1]-minmaxvec[maxind,1])
                        duNLmin = np.abs(SGuprime_O_u_tau[r,i-1]-minmaxvec[minind-1,0])
                        x_test_tensor = torch.tensor([[min(duNLmin, duNLmax)]], dtype=torch.float32)
                        pi, mu, sigma = Given_dist_NLmin_NLmax_model(x_test_tensor) 
                        randomvar = np.random.uniform(0,1)
                        duNN = 0
                        for j in range(Number_of_Gaussians):
                            duNN += pi[0,j]*norm.ppf(randomvar,loc = mu[0,j],scale = sigma[0,j])
                        if (SGuprime_O_u_tau[r,i-1]+np.exp(duNN)*res)>minmaxvec[maxind,1]:
                            SGuprime_O_u_tau[r,i] = minmaxvec[maxind,1]
                            maxind += 1
                            i += 1
                            is_min_sampled = not is_min_sampled
                        else:
                            SGuprime_O_u_tau[r,i] = SGuprime_O_u_tau[r,i-1]+np.exp(duNN)*res
                            i += 1
            sigma1 = std_uprime_O_u_tau_chopped[0,r]
            sigma2 = std_wprime_O_u_tau_chopped[0,r]
            SGuprime[r,:] = SGuprime_O_u_tau[r,:]*obj.u_tau
            SGuprime[r,:] = SGuprime - np.mean(SGuprime[r,:], axis = 1)
            SGwprime = rho_uw * (sigma2/sigma1) * SGuprime[r,:] + np.sqrt(1- rho_uw**2) * np.random.normal(0, sigma2, (1, SGuprime.shape[1]))
            SGwprime[r,:] = SGwprime[r,:] - np.mean(SGwprime[r,:], axis = 1)
        else:
            sigma1 = std_uprime_O_u_tau_chopped[0,r]
            sigma2 = std_wprime_O_u_tau_chopped[0,r]
            SGuprime[r,:] = 0.99*SGuprime[r-1,:] + np.sqrt(1-0.99**2)*np.random.normal(0, sigma1, (1, SGuprime.shape[1]))
            SGuprime[r,:] = SGuprime[r,:] - np.mean(SGuprime[r,:], axis = 1)

            SGwprime[r,:] = rho_uw * (sigma2/sigma1) * SGuprime[r,:] + np.sqrt(1- rho_uw**2) * np.random.normal(0, sigma2, (1, SGuprime.shape[1]))
            SGwprime[r,:] = SGwprime[r,:] - np.mean(SGwprime[r,:], axis = 1)

    for S in range(obj.u.shape[2]):
        for r in range(obj.z.shape[0]):
            obj.u[r,:,S] += SGuprime[r, S*obj.x.shape[0]:(S+1)*obj.x.shape[0]]
            obj.w[r,:,S] += SGwprime[r, S*obj.x.shape[0]:(S+1)*obj.x.shape[0]]
