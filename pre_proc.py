import h5py
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt

"Load Generated data(Velocity fields)"

print("Start reading the Data...\n")

WT7VF = sio.loadmat('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VLFIELD(WindTunnel(m1)).mat')
WT10VF = sio.loadmat('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VLFIELD(WindTunnel(m2)).mat')
ASLVF = sio.loadmat('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VLFIELD(ASL).mat')


uGenWT7 = np.array(WT7VF.get('ufieldGenWT7'))
uprimeGenWT7 = np.array(WT7VF.get('uprimefieldGenWT7'))
wGenWT7 = np.array(WT7VF.get('wfieldGenWT7'))

uGenWT10 = np.array(WT10VF.get('ufieldGenWT10'))
uprimeGenWT10 = np.array(WT10VF.get('uprimefieldGenWT10'))
wGenWT10 = np.array(WT10VF.get('wfieldGenWT10'))

uGenASL = np.array(ASLVF.get('ufieldGenABL'))
uprimeGenASL = np.array(ASLVF.get('uprimefieldGenABL'))
wGenASL = np.array(ASLVF.get('wfieldGenABL'))

del WT7VF, WT10VF, ASLVF


"Load Properties of Generated field"

WT7Prof = sio.loadmat('G:/My Drive/Research/DATABASE/Generated UMZ/UMZ(WindTunnel(m1)).mat')
WT10Prof = sio.loadmat('G:/My Drive/Research/DATABASE/Generated UMZ/UMZ(WindTunnel(m2)).mat')
ASLProf = sio.loadmat('G:/My Drive/Research/DATABASE/Generated UMZ/UMZ(ASL).mat')


zGenWT7 = np.array(WT7Prof.get('meanzWT7'))
zGenWT10 = np.array(WT10Prof.get('meanzWT10'))
zGenASL = np.array(ASLProf.get('meanzABL'))

del WT7Prof, WT10Prof, ASLProf

"Load Experimental dataset"

PIVWT7 = h5py.File(r'G:/My Drive/Research/DATABASE/Experimental/Wind Tunnel/PIV/Mesh/7[m s]/mesh_07ms_a_velocity.mat', 'r')
PIVWT10 = h5py.File(r'G:/My Drive/Research/DATABASE/Experimental/Wind Tunnel/PIV/Mesh/10[m s]/mesh_10ms_a_velocity.mat', 'r')
HotWireWT7 = h5py.File(r'G:/My Drive/Research/DATABASE/Experimental/Wind Tunnel/Hotwire/7[m s]/mesh_07ms_a_Xwire.mat','r')
HotWireWT10 = h5py.File(r'G:/My Drive/Research/DATABASE/Experimental/Wind Tunnel/Hotwire/10[m s]/mesh_10ms_a_Xwire.mat','r')
SLPIV = h5py.File(r'G:/My Drive/Research/DATABASE/Experimental/ASL/SLPIV/SLPIV_120Hz_Velocities.mat', 'r')
# for key in SLPIV.keys():
#     print(key)
# SLPIV = SLPIV['PIV_full']
# for key in SLPIV.keys():
#     print(key)
SonicASL = sio.loadmat('G:\My Drive\Research\DATABASE\Experimental\ASL\Sonic\SonicUTD_2022-02-22')
SonicASL = SonicASL['Sonic']



uPIVWT7 = np.array(PIVWT7.get('u'))
wPIVWT7 = np.array(PIVWT7.get('w'))
xPIVWT7 = np.array(PIVWT7.get('x'))
zPIVWT7 = np.array(PIVWT7.get('z'))

uPIVWT10 = np.array(PIVWT10.get('u'))
wPIVWT10 = np.array(PIVWT10.get('w'))
xPIVWT10 = np.array(PIVWT10.get('x'))
zPIVWT10 = np.array(PIVWT10.get('z'))

uHotWT7 = np.array(HotWireWT7.get('u'))
wHotWT7 = np.array(HotWireWT7.get('w'))
zHotWT7 = np.array(HotWireWT7.get('z'))
fsHotWT7 = 10000

uHotWT10 = np.transpose(np.array(HotWireWT10.get('u')))
wHotWT10 = np.transpose(np.array(HotWireWT10.get('w')))
zHotWT10 = np.transpose(np.array(HotWireWT10.get('z')))
fsHotWT10 = 10000


uSLPIV = SLPIV['u']
uSLPIV = np.transpose(np.array(uSLPIV[:,28:80,9:58]),(2, 1, 0))
wSLPIV = SLPIV['w']
wSLPIV = np.transpose(np.array(wSLPIV[:,28:80,9:58]),(2, 1, 0))
zSLPIV = SLPIV['Z']
zSLPIV = np.array(zSLPIV[0,9:58])
xSLPIV = SLPIV['X']
xSLPIV = np.array(xSLPIV[28:80,0])
fsSLPIV = 120

uSonicASL = SonicASL['Ux']
uSonicASL = uSonicASL[0,0]
wSonicASL = SonicASL['Uz']
wSonicASL = wSonicASL[0,0]
vSonicASL = SonicASL['Uy']
vSonicASL = vSonicASL[0,0]
fsSonicASL = SonicASL['Fsamp']
fsSonicASL = fsSonicASL[0,0]

U1 = np.mean(uSonicASL, axis = 0)
V1 = np.mean(vSonicASL, axis = 0)
alfa = np.arctan(V1/U1)
u1 = uSonicASL*np.cos(alfa) + vSonicASL*np.sin(alfa)
v1 = -uSonicASL*np.sin(alfa) + vSonicASL*np.cos(alfa)
w1 = wSonicASL

U2 = np.mean(u1, axis = 0)
V2 = np.mean(v1, axis = 0)
W2 = np.mean(w1, axis = 0)


beta = np.arctan(W2/U2)
uSonicASL = u1*np.cos(beta) + w1*np.sin(beta)
vSonicASL = v1
wSonicASL = -u1*np.sin(beta)+ w1*np.cos(beta)



del PIVWT7, PIVWT10, HotWireWT7, HotWireWT10, SLPIV, SonicASL



"Modifications on Generated Field"

uGenWT7 = uGenWT7[:,1000:]
uprimeGenWT7 = uprimeGenWT7[:,1000:]
wGenWT7 = wGenWT7[:,1000:]

uGenWT10 = uGenWT10[:,1000:]
uprimeGenWT10 = uprimeGenWT10[:,1000:]
wGenWT10 = wGenWT10[:,1000:]

uGenASL = uGenASL[:,1000:]
uprimeGenASL = uprimeGenASL[:,1000:]
wGenASL = wGenASL[:,1000:]

"Low-pass filter for ASL datasets (SLPIV, Sonic)"

Ns_SLPIV = uSLPIV.shape[2]
TurnT = 120 #Trunover time scale[s]
n = int(np.ceil(Ns_SLPIV/(TurnT*fsSLPIV))) #Number of frequencies needed to be considered for mean velocity reconstruction

PSu_SLPIV = np.fft.fft(uSLPIV, axis = 2)
PSU_SLPIV = np.zeros(PSu_SLPIV.shape, dtype = complex)
PSU_SLPIV[:,:,0:n+1] = PSu_SLPIV[:,:,0:n+1]
USLPIV = np.real(np.fft.ifft(PSU_SLPIV, axis = 2))

# plt.figure()
# plt.plot(uSLPIV[25,25,:])
# plt.plot(USLPIV[25,25,:])
# plt.show()

PSw_SLPIV = np.fft.fft(wSLPIV, axis = 2)
PSW_SLPIV = np.zeros(PSw_SLPIV.shape, dtype = complex)
PSW_SLPIV[:,:,0:n+1] = PSw_SLPIV[:,:,0:n+1]
WSLPIV = np.real(np.fft.ifft(PSW_SLPIV, axis = 2))

# plt.figure()
# plt.plot(wSLPIV[25,25,:])
# plt.plot(WSLPIV[25,25,:])
# plt.show()

Ns_SonicASL = uSonicASL.shape[0]
n = int(np.ceil(Ns_SonicASL/(TurnT*fsSonicASL[0,0])))

PSu_SonicASL = np.fft.fft(uSonicASL, axis = 0)
PSU_SonicASL = np.zeros(PSu_SonicASL.shape, dtype = complex)
PSU_SonicASL[0:n+1,0] = PSu_SonicASL[0:n+1,0]
USonicASL = np.real(np.fft.ifft(PSU_SonicASL, axis = 0))

# plt.figure()
# plt.plot(uSonicASL)
# plt.plot(USonicASL)
# plt.show()

deltaWT7=0.4078648
deltaWT10=0.3904
deltaASL=93

kappa=0.39
AFR=8.5

u_tauWT7=0.39
u_tauWT10=0.56
u_tauASL=0.4

z0WT7=0.00062
z0WT10=0.000631
z0ASL=0.0033

nuWT7=1.57*10**-5
nuWT10=1.57*10**-5
nuASL=1.26*10**-5

ksWT7=z0WT7/np.exp(-kappa*AFR)
ksWT10=z0WT10/np.exp(-kappa*AFR)
ksABL=z0ASL/np.exp(-kappa*AFR)

DelxGenWT7 = 1.09*10**-2
DelxGenWT10 = 1.03*10**-2
DelxGenASL = 3.68*10**-1

print("Data reading is completed.\n")