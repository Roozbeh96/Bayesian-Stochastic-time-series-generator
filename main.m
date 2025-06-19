

%% Data prepration
clc
clear
[data] = data_prep();

%% Make objects: Experimental PIV and SLPIV


VF_PIVWT7 = VFfeaturedVorX("PIVWT7","Exp",data.uPIVWT7,data.uprimePIVWT7,data.wPIVWT7,data.xPIVWT7,...
                            data.zPIVWT7,data.deltaWT7,data.u_tauWT7,NaN,NaN,...
                            data.z0WT7,data.nuWT7,NaN);
VF_PIVWT10 = VFfeaturedVorX("PIVWT10","Exp",data.uPIVWT10,data.uprimePIVWT10,data.wPIVWT10,data.xPIVWT10,...
                            data.zPIVWT10,data.deltaWT10,data.u_tauWT10,NaN,NaN,...
                            data.z0WT10,data.nuWT10,NaN);
VF_SLPIVASL = VFfeaturedVorX("SLPIVASL","Exp",data.uSLPIVASL,data.uprimeSLPIVASL,data.wSLPIVASL,data.xSLPIVASL,...
                            data.zSLPIVASL,data.deltaASL,data.u_tauASL,NaN,NaN,...
                            data.z0ASL,data.nuASL,data.fsSLPIVASL);
%% Make objects: Experimental H-W
VF_HotWT7 = VFfeaturedVorX("H-WWT7","Exp",data.uHotWT7,data.uprimeHotWT7,data.wHotWT7,NaN,...
                            data.zHotWT7,data.deltaWT7,data.u_tauWT7,NaN,NaN,...
                            data.z0WT7,data.nuWT7,data.fsHotWT7);          
VF_HotWT10 = VFfeaturedVorX("H-WWT10","Exp",data.uHotWT10,data.uprimeHotWT10,data.wHotWT10,NaN,...
                            data.zHotWT10,data.deltaWT10,data.u_tauWT10,NaN,NaN,...
                            data.z0WT10,data.nuWT10,data.fsHotWT10);   
VF_SonicASL = VFfeaturedVorX("SonicASL","Exp",data.uSonicASL,data.uprimeSonicASL,data.wSonicASL,NaN,...
                            2,data.deltaASL,data.u_tauASL,NaN,NaN,...
                            data.z0ASL,data.nuASL,data.fsSonicASL);   
%% Make objects: Generated 

VF_GenWT7 = VFfeaturedVorX("GenWT7","Gen",data.uGenWT7,data.uprimeGenWT7,data.wGenWT7,NaN,...
                            data.zGenWT7,data.deltaWT7,data.u_tauWT7,1,data.DelxGenWT7,...
                            data.z0WT7,data.nuWT7, NaN);
VF_GenWT10 = VFfeaturedVorX("GenWT10","Gen",data.uGenWT10,data.uprimeGenWT10,data.wGenWT10,NaN,...
                            data.zGenWT10,data.deltaWT10,data.u_tauWT10,1,data.DelxGenWT10,...
                            data.z0WT7,data.nuWT10, NaN);
VF_GenASL = VFfeaturedVorX("GenASL","Gen",data.uGenASL,data.uprimeGenASL,data.wGenASL,NaN,...
                            data.zGenASL,data.deltaASL,data.u_tauASL,1,data.DelxGenASL,...
                            data.z0ASL,data.nuASL, NaN);
% VF_GenASL = VFfeaturedVorX("GenASL","Gen",zeros(size(data.uGenASL)),zeros(size(data.uprimeGenASL)),zeros(size(data.wGenASL)),NaN,...
%                             data.zGenASL,data.deltaASL,data.u_tauASL,1,data.DelxGenASL,...
%                             data.z0ASL,data.nuASL, NaN);
%% clear data
clear('data')
%% Lambda_{ci} for experimental

VF_PIVWT7.Lambdaci()
VF_PIVWT10.Lambdaci()
VF_SLPIVASL.Lambdaci()
%% Spectral analysis for experimental
VF_HotWT7.Spectral_analysis()
VF_HotWT10.Spectral_analysis()
VF_SonicASL.Spectral_analysis()
%% Structure Function for experimental
VF_HotWT7.Structure_function()
VF_HotWT10.Structure_function()
VF_SonicASL.Structure_function()

%% Lambda_T for experimental
VF_HotWT7.Lambda_T_computation()
VF_HotWT10.Lambda_T_computation()
VF_SonicASL.Lambda_T_computation()
%% divide the field into LSMs

VF_GenWT7.Increasing_Res_nonlinear(1)%54
VF_GenWT10.Increasing_Res_nonlinear(1)%54
VF_GenASL.Increasing_Res_nonlinear(1)%54

save('temp.mat', 'VF_GenWT7', 'VF_GenWT10', 'VF_GenASL');
data_temp = load('temp.mat');

X = data_temp.VF_GenWT7;
Y = data_temp.VF_GenWT10;
Z = data_temp.VF_GenASL;
%% Detecting shear layers

VF_GenWT7.detecting_shearlayers()
% VF_GenWT10.detecting_shearlayers()
% VF_GenASL.detecting_shearlayers()
%% Increasing Resolution(use 'makima' for the interpolation)

VF_GenWT7.Increasing_Res_nonlinear(32)%50
VF_GenWT10.Increasing_Res_nonlinear(25)%50
VF_GenASL.Increasing_Res_nonlinear(20)%55

%% Increasing Resolution(use stochastic approach for the interpolation)

VF_GenWT7.Increasing_Res_sto(32)%25
VF_GenWT10.Increasing_Res_sto(25)%25
VF_GenASL.Increasing_Res_sto(20)%55
%% Lambda_{ci} for Generated

VF_GenWT7.Lambdaci()
VF_GenWT10.Lambdaci()
VF_GenASL.Lambdaci()

%% Vortex modeling 1(fully stochastic, place at lambda_ci~=0)
% close all
VF_GenWT7.VortXmodeling1()
VF_GenWT10.VortXmodeling1()
VF_GenASL.VortXmodeling1()
%% Vortex modeling 2(radius generated stochastic, place at shear layer)
% close all
VF_GenWT7.VortXmodeling2()


%% Gaussian-kernel(First_filter)

% VF_GenWT7.Gauss_kernel(1,1,[VF_GenWT7.Delx*0.5 0;...
%                 0 VF_GenWT7.Delx*0.5])
% VF_GenWT10.Gauss_kernel(1,1,[VF_GenWT10.Delx*0.5 0;...
%                 0 VF_GenWT10.Delx*0.5])
% VF_GenASL.Gauss_kernel(1,1,[VF_GenASL.Delx*0.5 0;...
%                 0 VF_GenASL.Delx*0.5])
%% Gaussian-kernel(Last_filter)

% VF_GenWT7.Gauss_kernel(1,0.2, [VF_GenWT7.Delx*1 VF_GenWT7.Delx*1*sqrt(tand(12.9));...
%                 VF_GenWT7.Delx*1*sqrt(tand(12.9)) tand(13)*VF_GenWT7.Delx*1])
VF_GenWT7.Gauss_kernel(0.4,0.4, [VF_GenWT7.Delx*0.00001 0;...
    0 VF_GenWT7.Delx*0.00001])
VF_GenWT10.Gauss_kernel(0.4,0.4, [VF_GenWT10.Delx*0.00001 0;...
    0 VF_GenWT10.Delx*0.00001])
VF_GenASL.Gauss_kernel(0.4,0.4, [VF_GenASL.Delx*0.00001 0;...
    0 VF_GenASL.Delx*0.00001])   
%% Structure Function
VF_GenWT7.Structure_function()
VF_GenWT10.Structure_function()
VF_GenASL.Structure_function()

%% Spectral analysis
VF_GenWT7.Spectral_analysis()
VF_GenWT10.Spectral_analysis()
VF_GenASL.Spectral_analysis()

%% Lambda_T
VF_GenWT7.Lambda_T_computation(VF_HotWT7)
VF_GenWT10.Lambda_T_computation(VF_HotWT10)
VF_GenASL.Lambda_T_computation(VF_SonicASL)


%% Two-point correlation
% VF_PIVWT7.Twop_corr(0.05)
VF_GenWT7.Twop_corr(0.05)


%% Saving
save('C:\Users\ehsan010\Desktop\VF_GenWT7(Beforepassingkernelres32makima-2).mat','VF_GenWT7')
save('C:\Users\ehsan010\Desktop\VF_GenWT10(Beforepassingkernelres25makima-2).mat','VF_GenWT10')
save('C:\Users\ehsan010\Desktop\VF_GenASL(Beforepassingkernelres20makima-2).mat','VF_GenASL')
VF_PIVWT7 = struct(VF_PIVWT7);
save('VF_PIVWT7.mat','VF_PIVWT7')
save('VF_PIVWT10.mat','VF_PIVWT10')
save('VF_SLPIVASL.mat','VF_SLPIVASL')
save('VF_HotWT7.mat','VF_HotWT7')
save('VF_HotWT10.mat','VF_HotWT10')
save('VF_SonicASL.mat','VF_SonicASL')