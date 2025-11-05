import h5py
import numpy as np
from Bayesian_ST.src.bayesian_st.Given_max_model import Given_max_model
from Given_min_model import MDN, Given_min_model
from scipy.stats import norm

def velocity_fluctuating_generation(obj):

    data1 = h5py.File(r'std_uprime_chopped.mat','r')
    data2 = h5py.File(r'std_wprime_chopped.mat','r')
    std_uprime_O_u_tau_chopped = np.array(data1.get('std_uprime_O_u_tau_chopped')).T
    std_uprime_O_u_tau_chopped = np.array(data2.get('std_wprime_O_u_tau_chopped')).T
    ''' 
    The correlation coeficient between u' and w' can be estimated using
    the chopped_signal_analysis.m line 80-84
    '''
    Number_of_Gaussians = 5
    rho_uw = -0.3 
    Given_min_model = MDN(1,1,Number_of_Gaussians)
    Given_min_model.load_state_dict(torch.load('G:\My Drive\Research\VFfeaturedVorX\Given_max_model_Long.pth'))
    Given_min_model.eval()
    x_test_tensor = torch.tensor([[x_test]], dtype=torch.float32)
    pi, mu, sigma = Given_min_model(x_test_tensor)

    Given_max_model = MDN(1,1,Number_of_Gaussians)
    Given_max_model.load_state_dict(torch.load('G:\My Drive\Research\VFfeaturedVorX\Given_max_model_Long.pth'))
    Given_max_model.eval()
    x_test_tensor = torch.tensor([[x_test]], dtype=torch.float32)
    pi, mu, sigma = Given_max_model(x_test_tensor)

    Given_dist_NLmin_NLmax_model = MDN(1,1,Number_of_Gaussians)
    Given_dist_NLmin_NLmax_model.load_state_dict(torch.load('G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model_Long.pth'))
    Given_dist_NLmin_NLmax_model.eval()
    x_test_tensor = torch.tensor([[x_test]], dtype=torch.float32)
    pi, mu, sigma = Given_dist_NLmin_NLmax_model(x_test_tensor) 

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
                    minmaxvec[minind,0] += pi[1,j]*norm.ppf(randomvar,loc = mu[0,j],scale = sigma[0,j])
            minind += 1
        is_min_sampled = not is_min_sampled
    
    minind = 0
    maxind = 0
    SGuprime = np.zeros((obj.z.shape[0], obj.x.shape[0], obj.u.shape[2]))
    SGwprime = np.zeros((obj.z.shape[0], obj.x.shape[0], obj.u.shape[2]))
    for S in range(obj.u.shape[2]):
        for r in range(obj.z.shape[0]):
            if r == 0:
                
                pass
            else:

                pass
