function [data] = data_prep()

fprintf('Reading data started...\n')
data = struct;
%<<<<Loading Generated modal velocity field>>>>
WT7VF = load('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VLFIELD(WindTunnel(m1)).mat');
WT10VF = load('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VLFIELD(WindTunnel(m2)).mat');
ASLVF = load('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VLFIELD(ASL).mat');

% WT7VF = load('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VFFIELD(WindTunnel(m1))_woSHlayer.mat');
% WT10VF = load('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VFFIELD(WindTunnel(m2))_woSHlayer.mat');
% ASLVF = load('G:\My Drive\Research\DATABASE\Generated, ordered,spaced VF\VFFIELD(ASL)_woSHlayer.mat');

WT7Prof = load('G:\My Drive\Research\DATABASE\Generated UMZ\UMZ(WindTunnel(m1)).mat');
WT10Prof = load('G:\My Drive\Research\DATABASE\Generated UMZ\UMZ(WindTunnel(m2)).mat');
ASLProf = load('G:\My Drive\Research\DATABASE\Generated UMZ\UMZ(ASL).mat');

uGenWT7 = WT7VF.ufieldGenWT7;
uprimeGenWT7 = uGenWT7-mean(uGenWT7,2);%WT7VF.uprimefieldGenWT7;%
wGenWT7 = WT7VF.wfieldGenWT7-mean(WT7VF.wfieldGenWT7,2);
zGenWT7 = WT7Prof.meanzWT7;

uGenWT10 = WT10VF.ufieldGenWT10;
uprimeGenWT10 = uGenWT10-mean(uGenWT10,2);%WT10VF.uprimefieldGenWT10;%
wGenWT10 = WT10VF.wfieldGenWT10-mean(WT10VF.wfieldGenWT10,2);
zGenWT10 = WT10Prof.meanzWT10;


uGenASL = ASLVF.ufieldGenABL;
uprimeGenASL = uGenASL-mean(uGenASL,2);%ASLVF.uprimefieldGenABL;%
wGenASL = ASLVF.wfieldGenABL-mean(ASLVF.wfieldGenABL,2);
zGenASL = ASLProf.meanzABL;

clear('WT7VF','WT10VF','ASLVF',...
    'WT7Prof', 'WT10Prof', 'ASLProf')


%<<<< Load experimental dataset >>>>


PIVWT7 = load('G:\My Drive\Research\DATABASE\Experimental\Wind Tunnel\PIV\Mesh\7[m s]\mesh_07ms_a_velocities.mat');
PIVWT10 = load('G:\My Drive\Research\DATABASE\Experimental\Wind Tunnel\PIV\Mesh\10[m s]\mesh_10ms_a_velocities.mat');
HotWireWT7 = load('G:\My Drive\Research\DATABASE\Experimental\Wind Tunnel\Hotwire\7[m s]\mesh_07ms_a_Xwire');
HotWireWT10 = load('G:\My Drive\Research\DATABASE\Experimental\Wind Tunnel\Hotwire\10[m s]\mesh_10ms_a_Xwire');
SLPIV = load('G:\My Drive\Research\DATABASE\Experimental\ASL\SLPIV\SLPIV_120Hz_Velocities.mat');
SonicASL = load('G:\My Drive\Research\DATABASE\Experimental\ASL\Sonic\SonicUTD_2022-02-22');


uPIVWT7 = PIVWT7.u;
uprimePIVWT7 = uPIVWT7-mean(mean(uPIVWT7,3),2);
wPIVWT7 = PIVWT7.w - mean(mean(PIVWT7.w,3),2);
xPIVWT7 = PIVWT7.x;
zPIVWT7 = PIVWT7.z;


uPIVWT10 = PIVWT10.u;
uprimePIVWT10 = uPIVWT10-mean(mean(uPIVWT10,3),2);
wPIVWT10 = PIVWT10.w - mean(mean(PIVWT10.w,3),2);
xPIVWT10 = PIVWT10.x;
zPIVWT10 = PIVWT10.z;


uHotWT7 = HotWireWT7.u;
uprimeHotWT7=uHotWT7-mean(uHotWT7,2);
wHotWT7 = HotWireWT7.w - mean(HotWireWT7.w,2);
wprimeHotWT7 = wHotWT7;
zHotWT7 = HotWireWT7.z;
fsHotWT7 = 10000;

uHotWT10 = HotWireWT10.u;
uprimeHotWT10 = uHotWT10-mean(uHotWT10,2);
wHotWT10 = HotWireWT10.w -mean(HotWireWT10.w,2);
wprimeHotWT10 = wHotWT10;
zHotWT10 = HotWireWT10.z;
fsHotWT10 = 10000;


uSLPIVASL=SLPIV.PIV_full.u(10:58,29:80,:);
wSLPIVASL = SLPIV.PIV_full.w(10:58,29:80,:);
zSLPIVASL=SLPIV.PIV_full.Z(10:58,1);
xSLPIVASL=SLPIV.PIV_full.X(1,29:80)';
fsSLPIVASL = 120;

uSonicASL = SonicASL.Sonic.Ux;
wSonicASL = SonicASL.Sonic.Uz;
vSonicASL = SonicASL.Sonic.Uy;
fsSonicASL = SonicASL.Sonic.Fsamp;

%<<<< Reprojecting Sonic dataset >>>>
U_m = mean(uSonicASL,1);
V_m = mean(vSonicASL,1);
alfa = atan(V_m./U_m);
u1 = uSonicASL*cos(alfa) + vSonicASL*sin(alfa);
v1 = -uSonicASL*sin(alfa) + vSonicASL*cos(alfa);
w1 = wSonicASL;

W1_m = mean(w1,1);
U1_m = mean(u1,1);

beta = atan(W1_m./U1_m);
uSonicASL = u1*cos(beta) + w1*sin(beta);
vSonicASL = v1;
wSonicASL = -u1*sin(beta)+ w1*cos(beta);

clear('PIVWT7','PIVWT10','HotWireWT7','HotWireWT10','SLPIV','SonicASL')


%<<<< Low-pass filter for Sonic and SLPIV >>>>

%SLPIV
Ns_SLPIVASL = size(uSLPIVASL,3);
TurnT = 120;
n = ceil(Ns_SLPIVASL/(TurnT*fsSLPIVASL));

PSu_SLPIVASL = fft(uSLPIVASL,[],3);
PSU_SLPIVASL = zeros(size(PSu_SLPIVASL));
PSU_SLPIVASL(:,:,1:n) = PSu_SLPIVASL(:,:,1:n);
USLPIVASL=real(ifft(PSU_SLPIVASL,[],3));


PSw_SLPIVASL = fft(wSLPIVASL,[],3);
PSW_SLPIVASL = zeros(size(PSw_SLPIVASL));
PSW_SLPIVASL(:,:,1:n) = PSw_SLPIVASL(:,:,1:n);
WSLPIVASL=real(ifft(PSW_SLPIVASL,[],3));


%Sonic

Ns_SonicASL = size(uSonicASL,1);
n = ceil(Ns_SonicASL/(TurnT*fsSonicASL));

PSu_SonicASL=fft(uSonicASL,[],1);
PSU_SonicASL = zeros(size(PSu_SonicASL));
PSU_SonicASL(1:n,1) = PSu_SonicASL(1:n,1);
USonicASL=real(ifft(PSU_SonicASL,[],1));

PSw_SonicASL=fft(wSonicASL,[],1);
PSW_SonicASL = zeros(size(PSw_SonicASL));
PSW_SonicASL(1:n,1) = PSw_SonicASL(1:n,1);
WSonicASL=real(ifft(PSW_SonicASL,[],1));

%<<<<---->>>>

uprimeSLPIVASL = uSLPIVASL-USLPIVASL;
wprimeSLPIVASL = wSLPIVASL-WSLPIVASL;

uprimeSonicASL = uSonicASL - USonicASL;
wprimeSonicASL = wSonicASL - WSonicASL;


%<<<< Removing first 1000 profiles >>>>
uGenWT7 = uGenWT7(:,1001:end);
uprimeGenWT7 = uprimeGenWT7(:,1001:end);
wGenWT7 = wGenWT7(:,1001:end);

uGenWT10 = uGenWT10(:,1001:end);
uprimeGenWT10 = uprimeGenWT10(:,1001:end);
wGenWT10 = wGenWT10(:,1001:end);

uGenASL = uGenASL(:,1001:end);
uprimeGenASL = uprimeGenASL(:,1001:end);
wGenASL = wGenASL(:,1001:end);




%<<<< Physical property of TBL >>>>

data.deltaWT7=0.4078648;
data.deltaWT10=0.3904;
data.deltaASL=93;

data.kappa=0.39;
data.AFR=8.5;

data.u_tauWT7=0.39;
data.u_tauWT10=0.56;
data.u_tauASL=0.4;

data.z0WT7=0.00062;
data.z0WT10=0.000631;
data.z0ASL=0.0033;

data.nuWT7=1.57*1e-5;
data.nuWT10=1.57*1e-5;
data.nuASL=1.26*1e-5;

data.ksWT7=data.z0WT7/exp(-data.kappa*data.AFR);
data.ksWT10=data.z0WT10/exp(-data.kappa*data.AFR);
data.ksABL=data.z0ASL/exp(-data.kappa*data.AFR);

data.DelxGenWT7 = 1.09*1e-2;
data.DelxGenWT10 = 1.03*1e-2;
data.DelxGenASL = 3.68*1e-1;

data.uGenWT7 = uGenWT7;
data.uprimeGenWT7 = uprimeGenWT7;
data.wGenWT7 = wGenWT7;
data.zGenWT7 = zGenWT7;
% data.Lambda_cirmsGenWT7 = 18;

data.uGenWT10 = uGenWT10;
data.uprimeGenWT10 = uprimeGenWT10;
data.wGenWT10 = wGenWT10;
data.zGenWT10 = zGenWT10;
% data.Lambda_cirmsGenWT10 = 0;

data.uGenASL = uGenASL;
data.uprimeGenASL = uprimeGenASL;
data.wGenASL = wGenASL;
data.zGenASL = zGenASL;
% data.Lambda_cirmsGenASL = 0;

data.uPIVWT7 = uPIVWT7;
data.uprimePIVWT7 = uprimePIVWT7;
data.wPIVWT7 = wPIVWT7;
data.xPIVWT7 = xPIVWT7;
data.zPIVWT7 = zPIVWT7;
% data.Lambda_cirmsPIVWT7 = 59;

data.uPIVWT10 = uPIVWT10;
data.uprimePIVWT10 = uprimePIVWT10;
data.wPIVWT10 = wPIVWT10;
data.xPIVWT10 = xPIVWT10;
data.zPIVWT10 = zPIVWT10;
% data.Lambda_cirmsPIVWT10 = 0;

data.uHotWT7 = uHotWT7;
data.uprimeHotWT7 = uprimeHotWT7;
data.wHotWT7 = wHotWT7;
data.wprimeHotWT7 = wprimeHotWT7;
data.zHotWT7 = zHotWT7;
data.fsHotWT7 = fsHotWT7;
data.etaHotWT7 = 0.000207199651618147;%%Based on D_{11}

data.uHotWT10 = uHotWT10;
data.uprimeHotWT10 = uprimeHotWT10;
data.wHotWT10 = wHotWT10;
data.wprimeHotWT10 = wprimeHotWT10;
data.zHotWT10 = zHotWT10;
data.fsHotWT10 = fsHotWT10;
data.etaHotWT10 = 0.000159751662919093;%%Based on D_{11}

data.uSLPIVASL = uSLPIVASL;
data.uprimeSLPIVASL = uprimeSLPIVASL;
data.wSLPIVASL = wprimeSLPIVASL;
data.wprimeSLPIVASL = wprimeSLPIVASL;
data.zSLPIVASL = zSLPIVASL;
data.xSLPIVASL = xSLPIVASL;
data.fsSLPIVASL = fsSLPIVASL;
data.etaSLPIVASL = 0.000455248515349699;%%Based on D_{11} and P=$\epsilon$
% data.Lambda_cirmsSLPIVASL = 0;

data.uSonicASL = uSonicASL;
data.uprimeSonicASL = uprimeSonicASL;
data.wSonicASL = wprimeSonicASL;
data.wprimeSonicASL = wprimeSonicASL;
data.vSonicASL = vSonicASL;
data.fsSonicASL = fsSonicASL;
data.etaSonicASL = 0.000390320092328079;%%Based on D_{11}


fprintf('Reading data finished\n')
end

