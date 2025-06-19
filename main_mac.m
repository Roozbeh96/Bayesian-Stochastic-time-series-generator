

%% Data prepration
clc
clear

[data] = data_prep_mac();

%% Make objects: Experimental(mac)


VF_PIVWT7 = VFfeaturedVorX_mac("PIVWT7","Exp",data.uPIVWT7,data.uprimePIVWT7,data.wPIVWT7,data.xPIVWT7,...
                            data.zPIVWT7,data.deltaWT7,data.u_tauWT7,NaN,NaN,data.Lambda_cirmsPIVWT7,...
                            data.z0WT7);
VF_PIVWT10 = VFfeaturedVorX_mac("PIVWT10","Exp",data.uPIVWT10,data.uprimePIVWT10,data.wPIVWT10,data.xPIVWT10,...
                            data.zPIVWT10,data.deltaWT10,data.u_tauWT10,NaN,NaN,data.Lambda_cirmsPIVWT10,...
                            data.z0WT10);
VF_SLPIVASL = VFfeaturedVorX_mac("SLPIVASL","Exp",data.uSLPIVASL,data.uprimeSLPIVASL,data.wSLPIVASL,data.xSLPIVASL,...
                            data.zSLPIVASL,data.deltaASL,data.u_tauASL,NaN,NaN,data.Lambda_cirmsSLPIVASL,...
                            data.z0ASL);

%% Make objects: Generated_mac

VF_GenWT7 = VFfeaturedVorX_mac("GenWT7","Gen",data.uGenWT7,data.uprimeGenWT7,data.wGenWT7,NaN,...
                            data.zGenWT7,data.deltaWT7,data.u_tauWT7,1,data.DelxGenWT7,data.Lambda_cirmsGenWT7,...
                            data.z0WT7);
VF_GenWT10 = VFfeaturedVorX_mac("GenWT10","Gen",data.uGenWT10,data.uprimeGenWT10,data.wGenWT10,NaN,...
                            data.zGenWT10,data.deltaWT10,data.u_tauWT10,1,data.DelxGenWT10,data.Lambda_cirmsGenWT10,...
                            data.z0WT7);
VF_GenASL = VFfeaturedVorX_mac("GenASL","Gen",data.uGenASL,data.uprimeGenASL,data.wGenASL,NaN,...
                            data.zGenASL,data.deltaASL,data.u_tauASL,1,data.DelxGenASL,data.Lambda_cirmsGenASL,...
                            data.z0WT7);

%% Lambda_{ci} for experimental

VF_PIVWT7.Lambdaci()
% VF_PIVWT10.Lambdaci()
% VF_SLPIVASL.Lambdaci()

%% Increasing Resolution

VF_GenWT7.Increasing_Res(54)
%% Gaussian-kernel

VF_GenWT7.Gauss_kernel(1,0.5)

%% Lambda_{ci} for Generated

VF_GenWT7.Lambdaci()
% VF_GenWT10.Lambdaci()
% VF_GenASL.Lambdaci()

%% Vortex modeling
% close all
VF_GenWT7.VortXmodeling()
%% Two-point correlation
% VF_PIVWT7.Twop_corr(0.05)
VF_GenWT7.Twop_corr(0.05)

%% Saving
save('VF_GenWT7(Beforepassingkernelmac).mat','VF_GenWT7')
VF_PIVWT7 = struct(VF_PIVWT7);
save('VF_PIVWT7(struc).mat','VF_PIVWT7')
