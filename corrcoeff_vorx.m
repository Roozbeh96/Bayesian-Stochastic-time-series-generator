%% Data preperation
Data1 = load('VF_PIVWT7.mat');
VF_PIVWT7 = Data1.VF_PIVWT7;

clear('Data1')

Data2 = load('VF_PIVWT10.mat');
VF_PIVWT10 = Data2.VF_PIVWT10;

clear('Data2')


Data3 = load('VF_SLPIVASL.mat');
VF_SLPIVASL = Data3.VF_SLPIVASL;

clear('Data3')

Data4 = load('VF_HotWT7.mat');
VF_HotWT7 = Data4.VF_HotWT7;
clear('Data4')

Data5 = load('VF_HotWT10.mat');
VF_HotWT10 = Data5.VF_HotWT10;
clear('Data5')

WT7VorX = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_07ms_vortex_properties.mat');
WT7VorX_scales = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_07ms_scales.mat');
corr_d_u_WT7 = corrcoef(WT7VorX.d_pro,-WT7VorX.dU_pro);

WT10VorX = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_10ms_vortex_properties.mat');
WT10VorX_scales = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_10ms_scales.mat');
corr_d_u_WT10 = corrcoef(WT10VorX.d_pro,-WT10VorX.dU_pro);

ASLVorX = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\asl_vortex_properties.mat');
ASLVorX_scales = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\asl_scales.mat');
corr_d_u_ASL = corrcoef(ASLVorX.d_pro,-ASLVorX.dU_pro);
%% WT7
corr_coef_local_WT7 = zeros();
ix = 1;

for r = 1:size(VF_PIVWT7.Lambda_ci, 1)
    VF_PIVWT7.Lambda_ci(r, :, :) = VF_PIVWT7.Lambda_ci(r, :, :) .* (abs(VF_PIVWT7.Lambda_ci(r, :, :)) >= 1.5*VF_PIVWT7.Lambda_cirms(r));
end

for S = 1: size(VF_PIVWT7.Lambda_ci,3)
Matrix = VF_PIVWT7.Lambda_ci(:,:,S);

binaryMatrix = Matrix ~= 0;

% Label connected components
[labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);  % 4-connectivity


% Calculate the mass center (centroid) for each connected component
for i = 1:numComponents
    % Extract the indices of the current component
    [rows, cols] = find(labeledMatrix == i);
    if size(rows,1)<30
        continue
    end
    values = Matrix(labeledMatrix == i);
    sign_region_mean = sign(mean(values));
    z_centroid = sum(VF_PIVWT7.z(rows) .* values) / sum(values);
    linear_indices = sub2ind(size(VF_PIVWT7.uprime), rows, cols, S*ones(size(rows,1),1));
    
    uprime = VF_PIVWT7.u(linear_indices)-mean(VF_PIVWT7.u(linear_indices),1);
    wprime = VF_PIVWT7.w(linear_indices)-mean(VF_PIVWT7.w(linear_indices),1);
    
%     uprime = VF_PIVWT7.uprime(linear_indices);
%     wprime = VF_PIVWT7.w(linear_indices);

    A = corrcoef(uprime,wprime);
    corr_coef_local_WT7(ix,1) = A(1,2);
    corr_coef_local_WT7(ix,2) = z_centroid;
    corr_coef_local_WT7(ix,3) = sign_region_mean;
    ix = ix + 1;
    
end
end
% histogram(corr_coef_local(:,1),200)
%% WT10
corr_coef_local_WT10 = zeros();
ix = 1;

for r = 1:size(VF_PIVWT10.Lambda_ci, 1)
    VF_PIVWT10.Lambda_ci(r, :, :) = VF_PIVWT10.Lambda_ci(r, :, :) .* (abs(VF_PIVWT10.Lambda_ci(r, :, :)) >= 1.5*VF_PIVWT10.Lambda_cirms(r));
end

for S = 1: size(VF_PIVWT10.Lambda_ci,3)
Matrix = VF_PIVWT10.Lambda_ci(:,:,S);

binaryMatrix = Matrix ~= 0;

% Label connected components
[labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);  % 4-connectivity


% Calculate the mass center (centroid) for each connected component
for i = 1:numComponents
    % Extract the indices of the current component
    [rows, cols] = find(labeledMatrix == i);
    if size(rows,1)<30
        continue
    end
    values = Matrix(labeledMatrix == i);
    sign_region_mean = sign(mean(values));
    z_centroid = sum(VF_PIVWT10.z(rows) .* values) / sum(values);
    linear_indices = sub2ind(size(VF_PIVWT10.uprime), rows, cols, S*ones(size(rows,1),1));
    
    uprime = VF_PIVWT10.u(linear_indices)-mean(VF_PIVWT10.u(linear_indices),1);
    wprime = VF_PIVWT10.w(linear_indices)-mean(VF_PIVWT10.w(linear_indices),1);
    
%     uprime = VF_PIVWT10.uprime(linear_indices);
%     wprime = VF_PIVWT10.w(linear_indices);

    A = corrcoef(uprime,wprime);
    corr_coef_local_WT10(ix,1) = A(1,2);
    corr_coef_local_WT10(ix,2) = z_centroid;
    corr_coef_local_WT10(ix,3) = sign_region_mean;
    ix = ix + 1;
    
end
end
%% ASL
corr_coef_local_ASL = zeros();
ix = 1;

for r = 1:size(VF_SLPIVASL.Lambda_ci, 1)
    VF_SLPIVASL.Lambda_ci(r, :, :) = VF_SLPIVASL.Lambda_ci(r, :, :) .* (abs(VF_SLPIVASL.Lambda_ci(r, :, :)) >= 1.5*VF_SLPIVASL.Lambda_cirms(r));
end

for S = 1: size(VF_SLPIVASL.Lambda_ci,3)
Matrix = VF_SLPIVASL.Lambda_ci(:,:,S);

binaryMatrix = Matrix ~= 0;

% Label connected components
[labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);  % 4-connectivity


% Calculate the mass center (centroid) for each connected component
for i = 1:numComponents
    % Extract the indices of the current component
    [rows, cols] = find(labeledMatrix == i);
    if size(rows,1)<30
        continue
    end
    values = Matrix(labeledMatrix == i);
    sign_region_mean = sign(mean(values));
    z_centroid = sum(VF_SLPIVASL.z(rows) .* values) / sum(values);
    linear_indices = sub2ind(size(VF_SLPIVASL.uprime), rows, cols, S*ones(size(rows,1),1));
    
    uprime = VF_SLPIVASL.u(linear_indices)-mean(VF_SLPIVASL.u(linear_indices),1);
    wprime = VF_SLPIVASL.w(linear_indices)-mean(VF_SLPIVASL.w(linear_indices),1);
    
%     uprime = VF_SLPIVASL.uprime(linear_indices);
%     wprime = VF_SLPIVASL.w(linear_indices);

    A = corrcoef(uprime,wprime);
    corr_coef_local_ASL(ix,1) = A(1,2);
    corr_coef_local_ASL(ix,2) = z_centroid;
    corr_coef_local_ASL(ix,3) = sign_region_mean;
    ix = ix + 1;
    
end
end
%% WT7,WT10,ASL
stat_rho_WT7 = zeros();
C1 = 10;
C2 = floor(size(VF_PIVWT7.z,1)/C1);
for i = 1:C2
    [r1] = find(corr_coef_local_WT7(:,2)>=VF_PIVWT7.z((i-1)*C1+1) &...
        corr_coef_local_WT7(:,2)<VF_PIVWT7.z((i)*C1) & corr_coef_local_WT7(:,3)<=0);
    tem_rho = corr_coef_local_WT7(r1,1);
    stat_rho_WT7(i,1) = mean(tem_rho,1);
    stat_rho_WT7(i,2) = std(tem_rho);
    stat_rho_WT7(i,3) = mean((VF_PIVWT7.z((i)*C1)+VF_PIVWT7.z((i-1)*C1+1))/2);
    [r2]= find(stat_rho_WT7(i,3)<=VF_HotWT7.z,1,'first');
    if r2 == 1
        stat_rho_WT7(i,4) = mean([VF_HotWT7.epsilon_str{r2,1}{:}],2);
    elseif isempty(r2)
        stat_rho_WT7(i,4)=NaN;
    else
        stat_rho_WT7(i,4) = ((VF_HotWT7.z(r2)-stat_rho_WT7(i,3))/(VF_HotWT7.z(r2)-VF_HotWT7.z(r2-1))*mean([VF_HotWT7.epsilon_str{r2,1}{:}],2)...
            +(stat_rho_WT7(i,3)-VF_HotWT7.z(r2-1))/(VF_HotWT7.z(r2)-VF_HotWT7.z(r2-1))*mean([VF_HotWT7.epsilon_str{r2-1,1}{:}],2));
    end
%     if i == size(VF_PIVWT7.z,1)/2
%         histogram(tem_rho,200)
%     end
end

stat_rho_WT10 = zeros();
C1 = 10;
C2 = floor(size(VF_PIVWT10.z,1)/C1);
for i = 1:C2
    [r1] = find(corr_coef_local_WT10(:,2)>=VF_PIVWT10.z((i-1)*C1+1) &...
        corr_coef_local_WT10(:,2)<VF_PIVWT10.z((i)*C1)& corr_coef_local_WT10(:,3)<=0);
    tem_rho = corr_coef_local_WT10(r1,1);
    stat_rho_WT10(i,1)= mean(tem_rho,1);
    stat_rho_WT10(i,2)=std(tem_rho);
    stat_rho_WT10(i,3)=mean((VF_PIVWT10.z((i)*C1)+VF_PIVWT10.z((i-1)*C1+1))/2);
    [r2]= find(stat_rho_WT10(i,3)<=VF_HotWT10.z,1,'first');
    if r2 == 1
        stat_rho_WT10(i,4) = mean([VF_HotWT10.epsilon_str{r2,1}{:}],2);
    elseif isempty(r2)
        stat_rho_WT10(i,4) = NaN;
    else
        stat_rho_WT10(i,4) = ((VF_HotWT10.z(r2)-stat_rho_WT10(i,3))/(VF_HotWT10.z(r2)-VF_HotWT10.z(r2-1))*mean([VF_HotWT10.epsilon_str{r2,1}{:}],2)...
            +(stat_rho_WT10(i,3)-VF_HotWT10.z(r2-1))/(VF_HotWT10.z(r2)-VF_HotWT10.z(r2-1))*mean([VF_HotWT10.epsilon_str{r2-1,1}{:}],2));
    end
%     if i == size(VF_PIVWT7.z,1)/2
%         histogram(tem_rho,200)
%     end
end

stat_rho_ASL = zeros();
C1 = 5;
C2 = floor(size(VF_SLPIVASL.z,1)/C1);
for i = 1:C2
    [r1] = find(corr_coef_local_ASL(:,2)>=VF_SLPIVASL.z((i-1)*C1+1) &...
        corr_coef_local_ASL(:,2)<VF_SLPIVASL.z((i)*C1));
    tem_rho = corr_coef_local_ASL(r1,1);
    stat_rho_ASL(i,1)= mean(tem_rho,1);
    stat_rho_ASL(i,2)=std(tem_rho);
    stat_rho_ASL(i,3)=mean((VF_SLPIVASL.z((i)*C1)+VF_SLPIVASL.z((i-1)*C1+1))/2);
%     if i == size(VF_PIVWT7.z,1)/2
%         histogram(tem_rho,200)
%     end
end
%%
figure
set(gcf,'Position',[565,444,474,356])
axes('Position',[0.09915611814346,0.160112359550562,0.873417721518986,0.806179775280899])


plot(stat_rho_WT7(:,1),stat_rho_WT7(:,3)/VF_PIVWT7.ks,'LineWidth',1.5,'color',...
    [0.17,1.00,0.72],'LineStyle','-','Marker','^','MarkerEdgeColor','k','MarkerFaceColor',...
    [0.17,1.00,0.72])
hold on
plot(stat_rho_WT10(:,1),stat_rho_WT10(:,3)/VF_PIVWT10.ks,'.','LineWidth',1.5,'color',...
    [0.00,0.30,0.60],'LineStyle','-','Marker','v','MarkerEdgeColor','k','MarkerFaceColor',...
    [0.00,0.30,0.60])
% plot(stat_rho_ASL(:,1),stat_rho_ASL(:,3)/VF_SLPIVASL.ks,'LineWidth',1.5,'color',...
%     'k','LineStyle','-')



plot([-0.15 -0.05],[2 2],'LineWidth',1.5,'color',...
    'r','LineStyle','--')
% plot([-0.15 -0.05],[2*VF_PIVWT10.ks/VF_PIVWT10.delta 2*VF_PIVWT10.ks/VF_PIVWT10.delta],'LineWidth',1.5,'color',...
%     'b','LineStyle','--')


set(gca,'TickLabelInterpreter','latex','FontSize',13,...
    'XGrid','on','YGrid','on','YScale','linear')
xlabel('$\mu_{\rho_{\mathrm{u,w}}^{\mathrm{VorX}}}$','Interpreter','Latex','FontSize',17,'position',...
    [-0.12499988079071,-0.676422752265178,-1]);
ylabel('z$/\mathrm{k}_{s}$','Interpreter','Latex','FontSize',17);
legend('PIV(m1)','PIV(m2)',... 
   'Interpreter','latex','FontSize',9,'Position',...
    [0.758660378031616,0.85325558131008,0.206399689317161,0.093539323967494],...
    'Numcolumns',1,'Orientation','Horizontal');
annotation(gcf,'textbox',...
    [0.298468354430376,0.311797752808989,0.155118143459917,0.092696629213483],...
    'String','z = $2\mathrm{k}_{s}$',...
    'Interpreter','latex',...
    'FontSize',13,...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([-0.25 0])
% ylim([0.0 8])

%%
figure
plot(stat_rho_WT7(:,2),stat_rho_WT7(:,3)/VF_PIVWT7.ks)
%% WT7
m = 2;

z = 1;
while stat_rho_WT7(z,3)/VF_PIVWT7.ks<=2
    if z==1
        Near_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2)>=0 &...
            corr_coef_local_WT7(:,2) <= stat_rho_WT7(z,3),:);
        Near_wall_d_pro = WT7VorX.d_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= stat_rho_WT7(z,3));
        Near_wall_U_omega_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= stat_rho_WT7(z,3));
        Near_wall_d_ret = WT7VorX.d_ret (WT7VorX.z_ret>=0 & WT7VorX.z_ret <= stat_rho_WT7(z,3));
        Near_wall_U_omega_ret_WT7 = WT7VorX.dU_ret (WT7VorX.z_ret>=0 & WT7VorX.z_ret <= stat_rho_WT7(z,3));
        
        Near_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro>=0 & WT7VorX.z_pro <= stat_rho_WT7(z,3))...
        - WT7VorX_scales.z.');
        [~, result_Near_wall_pro] = min(Near_wall_diffpro, [], 2);

        Near_wall_diffret = abs(WT7VorX.z_ret(WT7VorX.z_ret>=0 & WT7VorX.z_ret <= stat_rho_WT7(z,3))...
            - WT7VorX_scales.z.');
        [~, result_Near_wall_ret] = min(Near_wall_diffret, [], 2);

    else
        Near_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2)>=stat_rho_WT7(z-1,3) &...
            corr_coef_local_WT7(:,2) <= stat_rho_WT7(z,3),:);
        Near_wall_d_pro = WT7VorX.d_pro (WT7VorX.z_pro>=stat_rho_WT7(z-1,3) & WT7VorX.z_pro <= stat_rho_WT7(z,3));
        Near_wall_U_omega_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro>=stat_rho_WT7(z-1,3) & WT7VorX.z_pro <= stat_rho_WT7(z,3));
        Near_wall_d_ret = WT7VorX.d_ret (WT7VorX.z_ret>=stat_rho_WT7(z-1,3) & WT7VorX.z_ret <= stat_rho_WT7(z,3));
        Near_wall_U_omega_ret_WT7 = WT7VorX.dU_ret (WT7VorX.z_ret>=stat_rho_WT7(z-1,3) & WT7VorX.z_ret <= stat_rho_WT7(z,3));
        
        Near_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro>=stat_rho_WT7(z-1,3) & WT7VorX.z_pro <= stat_rho_WT7(z,3))...
        - WT7VorX_scales.z.');
        [~, result_Near_wall_pro] = min(Near_wall_diffpro, [], 2);

        Near_wall_diffret = abs(WT7VorX.z_ret(WT7VorX.z_ret>=stat_rho_WT7(z-1,3) & WT7VorX.z_ret <= stat_rho_WT7(z,3))...
            - WT7VorX_scales.z.');
        [~, result_Near_wall_ret] = min(Near_wall_diffret, [], 2);

    end
    Near_wall_vortices_WT7_pro = Near_wall_vortices_WT7(Near_wall_vortices_WT7(:,3)==-1,:);
    Near_wall_vortices_WT7_retro = Near_wall_vortices_WT7(Near_wall_vortices_WT7(:,3)==+1,:);
    
    [pdf_rho_Near_wall_vortices_pro,rho_Near_wall_vortices_pro] = ksdensity(Near_wall_vortices_WT7_pro(:,1),...
        linspace(min(Near_wall_vortices_WT7_pro(:,1)), max(Near_wall_vortices_WT7_pro(:,1)), 25),'NumPoints',25);
    [pdf_rho_Near_wall_vortices_retro,rho_Near_wall_vortices_retro] = ksdensity(Near_wall_vortices_WT7_retro(:,1),...
        linspace(min(Near_wall_vortices_WT7_retro(:,1)), max(Near_wall_vortices_WT7_retro(:,1)), 25),'NumPoints',25);

    [pdfr_proOlam_Near_wall,r_proOlamp_Near_wall] = ksdensity(0.5*Near_wall_d_pro./...
            WT7VorX_scales.lambda(result_Near_wall_pro),...
            'NumPoints',25);
    [pdfr_retOlam_Near_wall,r_retOlamp_Near_wall] = ksdensity(0.5*Near_wall_d_ret./...
            WT7VorX_scales.lambda(result_Near_wall_ret),...
            'NumPoints',25);
        
    [pdfU_omega_proOutau_Near_wall,dU_proOutaup_Near_wall] = ksdensity(-Near_wall_U_omega_pro_WT7/VF_PIVWT7.u_tau,...
        'NumPoints',25);
    [pdfU_omega_retOutau_Near_wall,dU_retOutaup_Near_wall] = ksdensity(Near_wall_U_omega_ret_WT7/VF_PIVWT7.u_tau,...
            'NumPoints',25);
        
    pdf_rho_Full_stat_Near_wall_vortices_WT7_pro(z,:) = pdf_rho_Near_wall_vortices_pro;
    rho_Full_stat_Near_wall_vortices_WT7_pro(z,:) = rho_Near_wall_vortices_pro;
    pdf_r_Full_stat_Near_wall_vortices_WT7_pro(z,:) = pdfr_proOlam_Near_wall;
    r_Full_stat_Near_wall_vortices_WT7_pro(z,:) = r_proOlamp_Near_wall;
    pdf_U_omega_Full_stat_Near_wall_vortices_WT7_pro(z,:) = pdfU_omega_proOutau_Near_wall;
    U_omega_Full_stat_Near_wall_vortices_WT7_pro(z,:) = dU_proOutaup_Near_wall;
    
    pdf_rho_Full_stat_Near_wall_vortices_WT7_retro(z,:) = pdf_rho_Near_wall_vortices_retro;
    rho_Full_stat_Near_wall_vortices_WT7_retro(z,:) = rho_Near_wall_vortices_retro;
    pdf_r_Full_stat_Near_wall_vortices_WT7_retro(z,:) = pdfr_retOlam_Near_wall;
    r_Full_stat_Near_wall_vortices_WT7_retro(z,:) = r_retOlamp_Near_wall;
    pdf_U_omega_Full_stat_Near_wall_vortices_WT7_retro(z,:) = pdfU_omega_retOutau_Near_wall;
    U_omega_Full_stat_Near_wall_vortices_WT7_retro(z,:) = dU_retOutaup_Near_wall;
    
    z = z+1;
end

Far_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2) > m*VF_PIVWT7.ks,:);
Far_wall_vortices_WT7_pro = Far_wall_vortices_WT7(Far_wall_vortices_WT7(:,3)==-1,:);
Far_wall_vortices_WT7_retro = Far_wall_vortices_WT7(Far_wall_vortices_WT7(:,3)==+1,:);
Farfrom_wall_d_pro = WT7VorX.d_pro (WT7VorX.z_pro > m*VF_PIVWT7.ks);
Farfrom_wall_dU_pro = WT7VorX.dU_pro (WT7VorX.z_pro > m*VF_PIVWT7.ks);
Farfrom_wall_d_ret = WT7VorX.d_ret (WT7VorX.z_ret > m*VF_PIVWT7.ks);
Farfrom_wall_dU_ret = WT7VorX.dU_ret (WT7VorX.z_ret > m*VF_PIVWT7.ks);

Farfrom_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro > m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Farfrom_wall_pro] = min(Farfrom_wall_diffpro, [], 2);

Farfrom_wall_diffret = abs(WT7VorX.z_ret(WT7VorX.z_ret > m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Farfrom_wall_ret] = min(Farfrom_wall_diffret, [], 2);


[pdfFar_wall_vortices_WT7_pro,rho_Far_wall_vortices_WT7_pro] = ksdensity(Far_wall_vortices_WT7_pro(:,1),...
    linspace(min(Far_wall_vortices_WT7_pro(:,1)), max(Far_wall_vortices_WT7_pro(:,1)), 25),'NumPoints',25);
[pdfFar_wall_vortices_WT7_retro,rho_Far_wall_vortices_WT7_retro] = ksdensity(Far_wall_vortices_WT7_retro(:,1),...
    linspace(min(Far_wall_vortices_WT7_retro(:,1)), max(Far_wall_vortices_WT7_retro(:,1)), 25),'NumPoints',25);

[pdfr_proOlam_Far_wall_WT7,r_proOlamp_Far_wall_WT7] = ksdensity(0.5*Farfrom_wall_d_pro./...
        WT7VorX_scales.lambda(result_Farfrom_wall_pro),...
        'NumPoints',25);
[pdfr_retOlam_Far_wall_WT7,r_retOlamp_Far_wall_WT7] = ksdensity(0.5*Farfrom_wall_d_ret./...
        WT7VorX_scales.lambda(result_Farfrom_wall_ret),...
        'NumPoints',25);

[pdfU_omega_proOutau_Far_wall_WT7,U_omega_proOutaup_Far_wall_WT7] = ksdensity(-Farfrom_wall_dU_pro/VF_PIVWT7.u_tau,...
        'NumPoints',25);
[pdfU_omega_retOutau_Far_wall_WT7,U_omega_retOutaup_Far_wall_WT7] = ksdensity(Farfrom_wall_dU_ret/VF_PIVWT7.u_tau,...
        'NumPoints',25);    


%% WT10

z = 1;
while stat_rho_WT10(z,3)/VF_PIVWT10.ks<=2
    if z==1
        Near_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2)>=0 &...
            corr_coef_local_WT10(:,2) <= stat_rho_WT10(z,3),:);
        Near_wall_d_pro = WT10VorX.d_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= stat_rho_WT10(z,3));
        Near_wall_U_omega_pro_WT7 = WT10VorX.dU_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= stat_rho_WT10(z,3));
        Near_wall_d_ret = WT10VorX.d_ret (WT10VorX.z_ret>=0 & WT10VorX.z_ret <= stat_rho_WT10(z,3));
        Near_wall_U_omega_ret_WT7 = WT10VorX.dU_ret (WT10VorX.z_ret>=0 & WT10VorX.z_ret <= stat_rho_WT10(z,3));
        
        Near_wall_diffpro = abs(WT10VorX.z_pro(WT10VorX.z_pro>=0 & WT10VorX.z_pro <= stat_rho_WT10(z,3))...
        - WT10VorX_scales.z.');
        [~, result_Near_wall_pro] = min(Near_wall_diffpro, [], 2);

        Near_wall_diffret = abs(WT10VorX.z_ret(WT10VorX.z_ret>=0 & WT10VorX.z_ret <= stat_rho_WT10(z,3))...
            - WT10VorX_scales.z.');
        [~, result_Near_wall_ret] = min(Near_wall_diffret, [], 2);
    else
        Near_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2)>=stat_rho_WT10(z-1,3) &...
            corr_coef_local_WT10(:,2) <= stat_rho_WT10(z,3),:);
        Near_wall_d_pro = WT10VorX.d_pro (WT10VorX.z_pro>=stat_rho_WT10(z-1,3) & WT10VorX.z_pro <= stat_rho_WT10(z,3));
        Near_wall_U_omega_pro_WT7 = WT10VorX.dU_pro (WT10VorX.z_pro>=stat_rho_WT10(z-1,3) & WT10VorX.z_pro <= stat_rho_WT10(z,3));
        Near_wall_d_ret = WT10VorX.d_ret (WT10VorX.z_ret>=stat_rho_WT10(z-1,3) & WT10VorX.z_ret <= stat_rho_WT10(z,3));
        Near_wall_U_omega_ret_WT7 = WT10VorX.dU_ret (WT10VorX.z_ret>=stat_rho_WT10(z-1,3) & WT10VorX.z_ret <= stat_rho_WT10(z,3));
        
        Near_wall_diffpro = abs(WT10VorX.z_pro(WT10VorX.z_pro>=stat_rho_WT10(z-1,3) & WT10VorX.z_pro <= stat_rho_WT10(z,3))...
        - WT10VorX_scales.z.');
        [~, result_Near_wall_pro] = min(Near_wall_diffpro, [], 2);

        Near_wall_diffret = abs(WT10VorX.z_ret(WT10VorX.z_ret>=stat_rho_WT10(z-1,3) & WT10VorX.z_ret <= stat_rho_WT10(z,3))...
            - WT10VorX_scales.z.');
        [~, result_Near_wall_ret] = min(Near_wall_diffret, [], 2);
    end
    Near_wall_vortices_WT10_pro = Near_wall_vortices_WT10(Near_wall_vortices_WT10(:,3)==-1,:);
    Near_wall_vortices_WT10_retro = Near_wall_vortices_WT10(Near_wall_vortices_WT10(:,3)==+1,:);
    [pdf_rho_Near_wall_vortices_pro,rho_Near_wall_vortices_pro] = ksdensity(Near_wall_vortices_WT10_pro(:,1),...
        linspace(min(Near_wall_vortices_WT10_pro(:,1)), max(Near_wall_vortices_WT10_pro(:,1)), 25),'NumPoints',25);
    [pdf_rho_Near_wall_vortices_retro,rho_Near_wall_vortices_retro] = ksdensity(Near_wall_vortices_WT10_retro(:,1),...
        linspace(min(Near_wall_vortices_WT10_retro(:,1)), max(Near_wall_vortices_WT10_retro(:,1)), 25),'NumPoints',25);
    
    [pdfr_proOlam_Near_wall,r_proOlamp_Near_wall] = ksdensity(0.5*Near_wall_d_pro./...
            WT10VorX_scales.lambda(result_Near_wall_pro),...
            'NumPoints',25);
    [pdfr_retOlam_Near_wall,r_retOlamp_Near_wall] = ksdensity(0.5*Near_wall_d_ret./...
            WT10VorX_scales.lambda(result_Near_wall_ret),...
            'NumPoints',25);
        
    [pdfU_omega_proOutau_Near_wall,dU_proOutaup_Near_wall] = ksdensity(-Near_wall_U_omega_pro_WT7/VF_PIVWT10.u_tau,...
        'NumPoints',25);
    [pdfU_omega_retOutau_Near_wall,dU_retOutaup_Near_wall] = ksdensity(Near_wall_U_omega_ret_WT7/VF_PIVWT10.u_tau,...
            'NumPoints',25);
    
    pdf_Full_stat_Near_wall_vortices_WT10_pro(z,:) = pdf_rho_Near_wall_vortices_pro;
    rho_Full_stat_Near_wall_vortices_WT10_pro(z,:) = rho_Near_wall_vortices_pro;
    pdf_r_Full_stat_Near_wall_vortices_WT10_pro(z,:) = pdfr_proOlam_Near_wall;
    r_Full_stat_Near_wall_vortices_WT10_pro(z,:) = r_proOlamp_Near_wall;
    pdf_U_omega_Full_stat_Near_wall_vortices_WT10_pro(z,:) = pdfU_omega_proOutau_Near_wall;
    U_omega_Full_stat_Near_wall_vortices_WT10_pro(z,:) = dU_proOutaup_Near_wall;
    
    pdf_Full_stat_Near_wall_vortices_WT10_retro(z,:) = pdf_rho_Near_wall_vortices_retro;
    rho_Full_stat_Near_wall_vortices_WT10_retro(z,:) = rho_Near_wall_vortices_retro;
    pdf_r_Full_stat_Near_wall_vortices_WT10_retro(z,:) = pdfr_retOlam_Near_wall;
    r_Full_stat_Near_wall_vortices_WT10_retro(z,:) = r_retOlamp_Near_wall;
    pdf_U_omega_Full_stat_Near_wall_vortices_WT10_retro(z,:) = pdfU_omega_retOutau_Near_wall;
    U_omega_Full_stat_Near_wall_vortices_WT10_retro(z,:) = dU_retOutaup_Near_wall;
    z = z+1;
end

Far_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2) > m*VF_PIVWT10.ks,:);
Far_wall_vortices_WT10_pro = Far_wall_vortices_WT10(Far_wall_vortices_WT10(:,3)==-1,:);
Far_wall_vortices_WT10_retro = Far_wall_vortices_WT10(Far_wall_vortices_WT10(:,3)==+1,:);
Farfrom_wall_d_pro = WT10VorX.d_pro (WT10VorX.z_pro > m*VF_PIVWT10.ks);
Farfrom_wall_dU_pro = WT10VorX.dU_pro (WT10VorX.z_pro > m*VF_PIVWT10.ks);
Farfrom_wall_d_ret = WT10VorX.d_ret (WT10VorX.z_ret > m*VF_PIVWT10.ks);
Farfrom_wall_dU_ret = WT10VorX.dU_ret (WT10VorX.z_ret > m*VF_PIVWT10.ks);

Farfrom_wall_diffpro = abs(WT10VorX.z_pro(WT10VorX.z_pro > m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Farfrom_wall_pro] = min(Farfrom_wall_diffpro, [], 2);

Farfrom_wall_diffret = abs(WT10VorX.z_ret(WT10VorX.z_ret > m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Farfrom_wall_ret] = min(Farfrom_wall_diffret, [], 2);

[pdfFar_wall_vortices_WT10_pro,rho_Far_wall_vortices_WT10_pro] = ksdensity(Far_wall_vortices_WT10_pro(:,1),...
    linspace(min(Far_wall_vortices_WT10_pro(:,1)), max(Far_wall_vortices_WT10_pro(:,1)), 25),'NumPoints',25);
[pdfFar_wall_vortices_WT10_retro,rho_Far_wall_vortices_WT10_retro] = ksdensity(Far_wall_vortices_WT10_retro(:,1),...
    linspace(min(Far_wall_vortices_WT10_retro(:,1)), max(Far_wall_vortices_WT10_retro(:,1)), 25),'NumPoints',25);

[pdfr_proOlam_Far_wall_WT10,r_proOlamp_Far_wall_WT10] = ksdensity(0.5*Farfrom_wall_d_pro./...
        WT10VorX_scales.lambda(result_Farfrom_wall_pro),...
        'NumPoints',25);
[pdfr_retOlam_Far_wall_WT10,r_retOlamp_Far_wall_WT10] = ksdensity(0.5*Farfrom_wall_d_ret./...
        WT10VorX_scales.lambda(result_Farfrom_wall_ret),...
        'NumPoints',25);

[pdfU_omega_proOutau_Far_wall_WT10,U_omega_proOutaup_Far_wall_WT10] = ksdensity(-Farfrom_wall_dU_pro/VF_PIVWT10.u_tau,...
        'NumPoints',25);
[pdfU_omega_retOutau_Far_wall_WT10,U_omega_retOutaup_Far_wall_WT10] = ksdensity(Farfrom_wall_dU_ret/VF_PIVWT10.u_tau,...
        'NumPoints',25); 

%% Fitting to both WT7 and WT10 for the radius(First near wall then Far)
m=2;

Near_wall_dU_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*VF_PIVWT7.ks);
Near_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*VF_PIVWT7.ks);
Near_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Near_wall_pro_WT7] = min(Near_wall_diffpro, [], 2);

Near_wall_dU_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*VF_PIVWT10.ks);
Near_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*VF_PIVWT10.ks);
Near_wall_diffpro = abs(WT10VorX.z_pro(WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Near_wall_pro_WT10] = min(Near_wall_diffpro, [], 2);

data_WT7 = double(0.5*Near_wall_d_pro_WT7./WT7VorX_scales.lambda(result_Near_wall_pro_WT7));
data_WT10 = double(0.5*Near_wall_d_pro_WT10./WT10VorX_scales.lambda(result_Near_wall_pro_WT10));




R1_WT7 = corrcoef(0.5*Near_wall_d_pro_WT7,-Near_wall_dU_pro_WT7);
R1_WT10 = corrcoef(0.5*Near_wall_d_pro_WT10,-Near_wall_dU_pro_WT10);
rho_near_du_r = abs(0.5*(R1_WT7(1,2)+R1_WT10(1,2))); % Target (Spearman) correlation %
n = 1e7; % Number of samples

randomvals = genGaussCop(rho_near_du_r, n);
rdelUw_near = randomvals(:,1);
rrw_near = randomvals(:,2);

% delUwOu_tau_near = -0.5*log(-rdelUw_near*(exp(-0.5/0.5)-exp(-4.5/0.5))+exp(-0.5/0.5));



Near_wall_U_omega_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*VF_PIVWT7.ks)/VF_PIVWT7.u_tau;
Near_wall_U_omega_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*VF_PIVWT10.ks)/VF_PIVWT10.u_tau;


data_WT7_near_pro = double(-Near_wall_U_omega_pro_WT7);
data_WT10_near_pro = double(-Near_wall_U_omega_pro_WT10);
a = 0.65;
b = 4.5;

data = [data_WT7_near_pro; data_WT10_near_pro];
data = data(data > a & data < b);

% Define Truncated Exponential PDF
exp_pdf_truncated = @(x, m, a, b) (1 ./ (m * (exp(-a / m) - exp(-b / m)))) .* exp(-x / m);

% Objective function (estimate only m)
objective_exp_truncated = @(m) -sum(log(exp_pdf_truncated(data, m, a, b)));

% Initial guess and bounds for m
initial_m = mean(data);
lb = 0.01;
ub = Inf;

% Optimize only for m
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
[opt_m, ~] = fmincon(objective_exp_truncated, initial_m, [], [], [], [], lb, ub, [], options);

% Compute C based on truncated normalization
C_opt = 1 / (opt_m * (exp(-a / opt_m) - exp(-b / opt_m)));

% Plot PDF over [0.5, 4.5]
x_values_near_u_omega = linspace(a, b, 100);
fitted_pdf_near_u_omega = exp_pdf_truncated(x_values_near_u_omega, opt_m , a, b);

delUwOu_tau_near = -opt_m*log(-rdelUw_near*(exp(-a/opt_m)-exp(-b/opt_m))+exp(-a/opt_m));

data = [data_WT7;data_WT10];
data = data(data>0.08);

% Log-Pareto PDF
log_pareto_pdf = @(x, mu, sigma, alpha, x_t) ...
    (x <= x_t) .* (1 ./ (x * sigma * sqrt(2 * pi)) .* exp(-((log(x) - mu).^2) ./ (2 * sigma^2))) + ...
    (x > x_t) .* (x_t^alpha ./ (x.^(alpha + 1)) .* ...
    (1 / (sigma * sqrt(2 * pi)) * exp(-((log(x_t) - mu)^2) / (2 * sigma^2))));

log_pareto_cdf = @(x, mu, sigma, alpha, x_t) ...
    (x <= x_t) .* normcdf(log(x), mu, sigma) + ...
    (x > x_t) .* (normcdf(log(x_t), mu, sigma) + ...
    (1 - normcdf(log(x_t), mu, sigma)) .* (1 - (x_t ./ x).^alpha));

% Objective function (negative log-likelihood)
objective_log_pareto = @(params) -sum(log(log_pareto_pdf(data, params(1), params(2), params(3), params(4))));

% Initial guess for parameters [mu, sigma, alpha, x_t]
initial_params = [mean(log(data)), std(log(data)),...
    5.5, 0.5];

% Bounds for the parameters [mu, sigma, alpha, x_t]
lb = [mean(log(data)), 0.01,5, prctile(data, 75)];%alpha = 5
ub = [0, Inf, 6, prctile(data, 90)];%alpha = 6

% Optimize using fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
[opt_params, ~] = fmincon(objective_log_pareto, initial_params, [], [], [], [], lb, ub, [], options);

% Extract optimized parameters
mu_opt = opt_params(1);
sigma_opt = opt_params(2);
alpha_opt = opt_params(3);
x_t_opt = opt_params(4);


x_values_near = logspace(log10(min(data)), log10(max(data)), 100);
fitted_pdf_near = log_pareto_pdf(x_values_near, mu_opt, sigma_opt,...
    alpha_opt,x_t_opt);
fitted_cdf_near = log_pareto_cdf(x_values_near, mu_opt, sigma_opt,...
    alpha_opt,x_t_opt);



F_xt = normcdf((log(x_t_opt) - mu_opt) / sigma_opt); % Log-normal CDF at x_t
rwOlambda_T_near = zeros(size(rrw_near));

%genearted close wall r_w/lambda For For u <= F_xt (Log-normal region)
idx_log_normal = rrw_near <= F_xt;
rwOlambda_T_near(idx_log_normal) = exp(mu_opt + sigma_opt * norminv(rrw_near(idx_log_normal)));

%genearted close wall r_w/lambda For For u > F_xt (Power-law region)
idx_power_law = rrw_near > F_xt;
rwOlambda_T_near(idx_power_law) = x_t_opt * (1 - (rrw_near(idx_power_law) - F_xt) / (1 - F_xt)).^(-1 / alpha_opt);




% Far

Far_wall_dU_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro > m*VF_PIVWT7.ks);
Far_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro > m*VF_PIVWT7.ks);
Far_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro > m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Far_wall_pro_WT7] = min(Far_wall_diffpro, [], 2);
            
Far_wall_dU_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro > m*VF_PIVWT10.ks);
Far_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro > m*VF_PIVWT10.ks);
Far_wall_diffpro = abs(WT10VorX.z_pro(WT10VorX.z_pro > m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Far_wall_pro_WT10] = min(Far_wall_diffpro, [], 2);

data_WT7 = double(0.5*Far_wall_d_pro_WT7./WT7VorX_scales.lambda(result_Far_wall_pro_WT7));
data_WT10 = double(0.5*Far_wall_d_pro_WT10./WT10VorX_scales.lambda(result_Far_wall_pro_WT10));




R2_WT7 = corrcoef(0.5*Far_wall_d_pro_WT7,-Far_wall_dU_pro_WT7);
R2_WT10 = corrcoef(0.5*Far_wall_d_pro_WT10,-Far_wall_dU_pro_WT10);
rho_Far_du_r = abs(0.5*(R2_WT7(1,2)+R2_WT10(1,2))); % Target (Spearman) correlation %
randomvals = genGaussCop(rho_Far_du_r, n);
rdelUw_Far = randomvals(:,1);
rrw_Far = randomvals(:,2);

% delUwOu_tau_Far = -0.37*log(-rdelUw_Far*(exp(-0.5/0.37)-exp(-4.5/0.37))+exp(-0.5/0.37));

Far_wall_U_omega_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro > m*VF_PIVWT7.ks)/VF_PIVWT7.u_tau;
Far_wall_U_omega_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro > m*VF_PIVWT10.ks)/VF_PIVWT10.u_tau;

data_WT7_Far_pro = double(-Far_wall_U_omega_pro_WT7);
data_WT10_Far_pro = double(-Far_wall_U_omega_pro_WT10);

a = 0.4;
b = 4.5;
data = [data_WT7_Far_pro ; data_WT10_Far_pro];
data = data(data > a & data < b);

% Define Truncated Exponential PDF
exp_pdf_truncated = @(x, m, a, b) (1 ./ (m * (exp(-a / m) - exp(-b / m)))) .* exp(-x / m);

% Objective function (estimate only m)
objective_exp_truncated = @(m) -sum(log(exp_pdf_truncated(data, m, a, b)));

% Initial guess and bounds for m
initial_m = mean(data);
lb = 0.01;
ub = Inf;

% Optimize only for m
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
[opt_m, ~] = fmincon(objective_exp_truncated, initial_m, [], [], [], [], lb, ub, [], options);

% Compute C based on truncated normalization
C_opt = 1 / (opt_m * (exp(-a / opt_m) - exp(-b / opt_m)));

% Plot PDF over [0.5, 4.5]
x_values_Far_u_omega = linspace(a, b, 100);
fitted_pdf_Far_u_omega = exp_pdf_truncated(x_values_Far_u_omega, opt_m, a, b);

delUwOu_tau_Far = -opt_m*log(-rdelUw_Far*(exp(-a/opt_m)-exp(-b/opt_m))+exp(-a/opt_m));

data = [data_WT7;data_WT10];
data = data(data>0.06);

% Log-Pareto PDF
log_pareto_pdf = @(x, mu, sigma, alpha, x_t) ...
    (x <= x_t) .* (1 ./ (x * sigma * sqrt(2 * pi)) .* exp(-((log(x) - mu).^2) ./ (2 * sigma^2))) + ...
    (x > x_t) .* (x_t^alpha ./ (x.^(alpha + 1)) .* ...
    (1 / (sigma * sqrt(2 * pi)) * exp(-((log(x_t) - mu)^2) / (2 * sigma^2))));

log_pareto_cdf = @(x, mu, sigma, alpha, x_t) ...
    (x <= x_t) .* normcdf(log(x), mu, sigma) + ...
    (x > x_t) .* (normcdf(log(x_t), mu, sigma) + ...
    (1 - normcdf(log(x_t), mu, sigma)) .* (1 - (x_t ./ x).^alpha));

% Objective function (negative log-likelihood)
objective_log_pareto = @(params) -sum(log(log_pareto_pdf(data, params(1), params(2), params(3), params(4))));

% Initial guess for parameters [mu, sigma, alpha, x_t]
initial_params = [mean(log(data)), std(log(data)),...
    5.5, 0.5];

% Bounds for the parameters [mu, sigma, alpha, x_t]
lb = [mean(log(data)), 0.01, 4.5, prctile(data, 75)];%alpha =4.5
ub = [0, Inf, 6, prctile(data, 90)];%alpha=6

% Optimize using fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
[opt_params, ~] = fmincon(objective_log_pareto, initial_params, [], [], [], [], lb, ub, [], options);

% Extract optimized parameters
mu_opt = opt_params(1);
sigma_opt = opt_params(2);
alpha_opt = opt_params(3);
x_t_opt = opt_params(4);

x_values_Far = logspace(log10(min(data)), log10(max(data)), 100);
fitted_pdf_Far = log_pareto_pdf(x_values_Far, mu_opt, sigma_opt,...
    alpha_opt,x_t_opt);
fitted_cdf_Far = log_pareto_cdf(x_values_Far, mu_opt, sigma_opt,...
    alpha_opt,x_t_opt);



F_xt = normcdf((log(x_t_opt) - mu_opt) / sigma_opt); % Log-normal CDF at x_t
rwOlambda_T_Far = zeros(size(rrw_Far));

% genearted far wall r_w/lambda For u <= F_xt (Log-normal region)
idx_log_normal = rrw_Far <= F_xt;
rwOlambda_T_Far(idx_log_normal) = exp(mu_opt + sigma_opt * norminv(rrw_Far(idx_log_normal)));

% genearted far wall r_w/lambda For u > F_xt (Power-law region)
idx_power_law = rrw_Far > F_xt;
rwOlambda_T_Far(idx_power_law) = x_t_opt * (1 - (rrw_Far(idx_power_law) - F_xt) / (1 - F_xt)).^(-1 / alpha_opt);



%% Plotting WT7
figure
for z = 1:size(pdf_rho_Full_stat_Near_wall_vortices_WT7_pro,1)
    plot(rho_Full_stat_Near_wall_vortices_WT7_pro(z,:),pdf_rho_Full_stat_Near_wall_vortices_WT7_pro(z,:),...
        'LineStyle','none','color','b','Marker','^');
        % Append the current legend entry to the list
    legend_entries{z} = sprintf('PIV(m1) $\\mathrm{z}/\\mathrm{k}_{s}$ = %.2f', stat_rho_WT7(z, 3) / VF_PIVWT7.ks);

    hold on
end

plot(rho_Far_wall_vortices_WT7_pro,pdfFar_wall_vortices_WT7_pro,...
        'LineStyle','none','color','r','Marker','^');
set(gca,'TickLabelInterpreter','latex','FontSize',13)
legend_entries{z+1} = 'PIV(m1) $\mathrm{z}/\mathrm{k}_{s} \geq 2$';
legend(legend_entries, 'Interpreter', 'latex');
xlabel('$\rho_{\mathrm{u,w}}^{\mathrm{VorX}}$','Interpreter','Latex','FontSize',17);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);

% figure
% for z = 1:size(pdf_rho_Full_stat_Near_wall_vortices_WT7_retro,1)
%     plot(rho_Full_stat_Near_wall_vortices_WT7_retro(z,:),pdf_rho_Full_stat_Near_wall_vortices_WT7_retro(z,:),...
%         'LineStyle','none','color','b','Marker','^');
%     hold on
% end
% plot(rho_Far_wall_vortices_WT7_retro,pdfFar_wall_vortices_WT7_retro,...
%         'LineStyle','none','color','r','Marker','^');


figure
for z = 1:size(pdf_r_Full_stat_Near_wall_vortices_WT7_pro,1)
    plot(r_Full_stat_Near_wall_vortices_WT7_pro(z,:),pdf_r_Full_stat_Near_wall_vortices_WT7_pro(z,:),...
        'LineStyle','none','color','b','Marker','^');
    legend_entries{z} = sprintf('PIV(m1) $\\mathrm{z}/\\mathrm{k}_{s}$ = %.2f', stat_rho_WT7(z, 3) / VF_PIVWT7.ks);
    hold on
end
plot(r_proOlamp_Far_wall_WT7,pdfr_proOlam_Far_wall_WT7,...
        'LineStyle','none','color','r','Marker','^');
plot(x_values_near,fitted_pdf_near,'k','LineWidth',1.5)
plot(x_values_Far,fitted_pdf_Far,'r','LineWidth',1.5)
set(gca,'TickLabelInterpreter','latex','FontSize',13,'YScale','log','XScale','log',...
    'XGrid','on','YGrid','on')
legend_entries{z+1} = 'PIV(m1) $\mathrm{z}/\mathrm{k}_{s} \geq 2$';
legend(legend_entries, 'Interpreter', 'latex');
xlabel('r$_{\omega}/\lambda_{T}$','Interpreter','Latex','FontSize',17);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
ylim([1e-5 5e1])
xlim([0.03 0.2e1])

% figure
% for z = 1:size(pdf_r_Full_stat_Near_wall_vortices_WT7_retro,1)
%     plot(r_Full_stat_Near_wall_vortices_WT7_retro(z,:),pdf_r_Full_stat_Near_wall_vortices_WT7_retro(z,:),...
%         'LineStyle','none','color','b','Marker','^');
%     hold on
% end
% plot(r_retOlamp_Far_wall_WT7,pdfr_retOlam_Far_wall_WT7,...
%         'LineStyle','none','color','r','Marker','^');
%   

figure
for z = 1:size(pdf_U_omega_Full_stat_Near_wall_vortices_WT7_pro,1)
    plot(U_omega_Full_stat_Near_wall_vortices_WT7_pro(z,:),pdf_U_omega_Full_stat_Near_wall_vortices_WT7_pro(z,:),...
        'LineStyle','none','color','b','Marker','^');
    hold on
end
[N,edges]=histcounts(data_WT7_near_pro,100, 'Normalization', 'pdf');

plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N ,...
        'LineStyle','none','color','k','Marker','O')
plot(U_omega_proOutaup_Far_wall_WT7,pdfU_omega_proOutau_Far_wall_WT7,...
        'LineStyle','none','color','r','Marker','^');
semilogy(x_values_near_u_omega,fitted_pdf_near_u_omega,'k','LineWidth',1.5)
semilogy(x_values_Far_u_omega,fitted_pdf_Far_u_omega,'r','LineWidth',1.5)
% semilogy((0.5:0.1:4.5),5.43*exp(-(0.5:0.1:4.5)./0.5),'k','LineWidth',1.5)
% semilogy((0.5:0.1:4.5),10.43*exp(-(0.5:0.1:4.5)./0.37),'r','LineWidth',1.5)
set(gca,'TickLabelInterpreter','latex','FontSize',13,'YScale','log',...
    'XGrid','on','YGrid','on')
legend_entries{z+1} = 'PIV(m1) $\mathrm{z}/\mathrm{k}_{s} \geq 2$';
legend(legend_entries, 'Interpreter', 'latex');
xlabel('u$_{\omega}/u_{\tau}$','Interpreter','Latex','FontSize',17);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
ylim([1e-7 1e1])
xlim([-0.5 5.5])    

% figure
% for z = 1:size(pdf_U_omega_Full_stat_Near_wall_vortices_WT7_retro,1)
%     plot(U_omega_Full_stat_Near_wall_vortices_WT7_retro(z,:),pdf_U_omega_Full_stat_Near_wall_vortices_WT7_retro(z,:),...
%         'LineStyle','none','color','b','Marker','^');
%     hold on
% end
% plot(U_omega_retOutaup_Far_wall_WT7,pdfU_omega_retOutau_Far_wall_WT7,...
%         'LineStyle','none','color','r','Marker','^');
    
%% Plotting WT10
figure
for z = 1:size(pdf_Full_stat_Near_wall_vortices_WT10_pro,1)
    plot(rho_Full_stat_Near_wall_vortices_WT10_pro(z,:),pdf_Full_stat_Near_wall_vortices_WT10_pro(z,:),...
        'LineStyle','none','color','b','Marker','^');
    hold on
    legend_entries{z} = sprintf('PIV(m2) $\\mathrm{z}/\\mathrm{k}_{s}$ = %.3f', stat_rho_WT10(z, 3) / VF_PIVWT10.ks);
end
plot(rho_Far_wall_vortices_WT10_pro,pdfFar_wall_vortices_WT10_pro,...
        'LineStyle','none','color','r','Marker','^'); 
set(gca,'TickLabelInterpreter','latex','FontSize',13)
legend_entries{z+1} = 'PIV(m2) $\mathrm{z}/\mathrm{k}_{s} \geq 2$';
legend(legend_entries, 'Interpreter', 'latex');
xlabel('$\rho_{\mathrm{u,w}}^{\mathrm{VorX}}$','Interpreter','Latex','FontSize',17);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);

% figure
% for z = 1:size(pdf_Full_stat_Near_wall_vortices_WT10_retro,1)
%     plot(rho_Full_stat_Near_wall_vortices_WT10_retro(z,:),pdf_Full_stat_Near_wall_vortices_WT10_retro(z,:),...
%         'LineStyle','none','color','b','Marker','^');
%     hold on
% end
% plot(rho_Far_wall_vortices_WT10_retro,pdfFar_wall_vortices_WT10_retro,...
%         'LineStyle','none','color','r','Marker','^');

figure
for z = 1:size(pdf_r_Full_stat_Near_wall_vortices_WT10_pro,1)
    plot(r_Full_stat_Near_wall_vortices_WT10_pro(z,:),pdf_r_Full_stat_Near_wall_vortices_WT10_pro(z,:),...
        'LineStyle','none','color','b','Marker','^');
    hold on
    legend_entries{z} = sprintf('PIV(m2) $\\mathrm{z}/\\mathrm{k}_{s}$ = %.3f', stat_rho_WT10(z, 3) / VF_PIVWT10.ks);
end
plot(r_proOlamp_Far_wall_WT10,pdfr_proOlam_Far_wall_WT10,...
        'LineStyle','none','color','r','Marker','^');
plot(x_values_near,fitted_pdf_near,'k','LineWidth',1.5)
plot(x_values_Far,fitted_pdf_Far,'r','LineWidth',1.5)       
set(gca,'TickLabelInterpreter','latex','FontSize',13,'YScale','log','XScale','log',...
    'XGrid','on','YGrid','on')
legend_entries{z+1} = 'PIV(m2) $\mathrm{z}/\mathrm{k}_{s} \geq 2$';
legend(legend_entries, 'Interpreter', 'latex');
xlabel('r$_{\omega}/\lambda_{T}$','Interpreter','Latex','FontSize',17);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
ylim([1e-5 5e1])
xlim([0.03 0.2e1])

% figure
% for z = 1:size(pdf_r_Full_stat_Near_wall_vortices_WT10_retro,1)
%     plot(r_Full_stat_Near_wall_vortices_WT10_retro(z,:),pdf_r_Full_stat_Near_wall_vortices_WT10_retro(z,:),...
%         'LineStyle','none','color','b','Marker','^');
%     hold on
% end
% plot(r_retOlamp_Far_wall_WT10,pdfr_retOlam_Far_wall_WT10,...
%         'LineStyle','none','color','r','Marker','^');
    
figure
for z = 1:size(pdf_U_omega_Full_stat_Near_wall_vortices_WT10_pro,1)
    plot(U_omega_Full_stat_Near_wall_vortices_WT10_pro(z,:),pdf_U_omega_Full_stat_Near_wall_vortices_WT10_pro(z,:),...
        'LineStyle','none','color','b','Marker','^');
    hold on
end
[N,edges]=histcounts(data_WT10_near_pro,100, 'Normalization', 'pdf');

plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N ,...
        'LineStyle','none','color','k','Marker','O')
plot(U_omega_proOutaup_Far_wall_WT10,pdfU_omega_proOutau_Far_wall_WT10,...
        'LineStyle','none','color','r','Marker','^');
semilogy(x_values_near_u_omega,fitted_pdf_near_u_omega,'k','LineWidth',1.5)
semilogy(x_values_Far_u_omega,fitted_pdf_Far_u_omega,'r','LineWidth',1.5)
% semilogy((0.5:0.1:4.5),5.43*exp(-(0.5:0.1:4.5)./0.5),'k','LineWidth',1.5)
% semilogy((0.5:0.1:4.5),10.43*exp(-(0.5:0.1:4.5)./0.37),'r','LineWidth',1.5)
set(gca,'TickLabelInterpreter','latex','FontSize',13,'YScale','log',...
    'XGrid','on','YGrid','on')
legend_entries{z+1} = 'PIV(m2) $\mathrm{z}/\mathrm{k}_{s} \geq 2$';
legend(legend_entries, 'Interpreter', 'latex');
xlabel('u$_{\omega}/u_{\tau}$','Interpreter','Latex','FontSize',17);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
ylim([1e-7 1e1])
xlim([-0.5 5.5]) 
% figure
% for z = 1:size(pdf_U_omega_Full_stat_Near_wall_vortices_WT10_retro,1)
%     plot(U_omega_Full_stat_Near_wall_vortices_WT10_retro(z,:),pdf_U_omega_Full_stat_Near_wall_vortices_WT10_retro(z,:),...
%         'LineStyle','none','color','b','Marker','^');
%     hold on
% end
% plot(U_omega_retOutaup_Far_wall_WT10,pdfU_omega_retOutau_Far_wall_WT10,...
%         'LineStyle','none','color','r','Marker','^');    
%%

figure;
histogram(delUwOu_tau_near, 'Normalization', 'pdf'); 
hold on;
[f, xi] = ksdensity(delUwOu_tau_near); 
plot(xi, f, 'r-', 'LineWidth', 2); 
set(gca,'YScale','log')

figure;
histogram(delUwOu_tau_Far, 'Normalization', 'pdf'); 
hold on;
[f, xi] = ksdensity(delUwOu_tau_Far); 
plot(xi, f, 'r-', 'LineWidth', 2); 
set(gca,'YScale','log')

figure;
histogram(rwOlambda_T_near, 'Normalization', 'pdf'); 
hold on;
[f, xi] = ksdensity(rwOlambda_T_near); 
plot(xi, f, 'r-', 'LineWidth', 2); 
set(gca,'YScale','log','XScale','log')
ylim([10^-5 10^1])

figure;
histogram(rwOlambda_T_Far, 'Normalization', 'pdf'); 
hold on;
[f, xi] = ksdensity(rwOlambda_T_Far); 
plot(xi, f, 'r-', 'LineWidth', 2); 
set(gca,'YScale','log','XScale','log')
ylim([10^-5 10^1])
%% Save the generated random numbers

save('delUwOu_tau_near.mat','delUwOu_tau_near')
save('delUwOu_tau_Far.mat','delUwOu_tau_Far')
save('rwOlambda_T_near.mat','rwOlambda_T_near')
save('rwOlambda_T_Far.mat','rwOlambda_T_Far')

load('delUwOu_tau_near.mat')
load('delUwOu_tau_Far.mat')
load('rwOlambda_T_near.mat')
load('rwOlambda_T_Far.mat')
%% ASL
z = 1;
while stat_rho_ASL(z,3)/VF_SLPIVASL.ks<=2
    if z==1
        Near_wall_vortices_ASL = corr_coef_local_ASL (corr_coef_local_ASL(:,2)>=0 &...
            corr_coef_local_ASL(:,2) <= stat_rho_ASL(z,3),:);
    else
        Near_wall_vortices_ASL = corr_coef_local_ASL (corr_coef_local_ASL(:,2)>=stat_rho_ASL(z-1,3) &...
            corr_coef_local_ASL(:,2) <= stat_rho_ASL(z,3),:);
    end
    Near_wall_vortices_ASL_pro = Near_wall_vortices_ASL(Near_wall_vortices_ASL(:,3)==-1,:);
    Near_wall_vortices_ASL_retro = Near_wall_vortices_ASL(Near_wall_vortices_ASL(:,3)==+1,:);
    [pdf_rho_Near_wall_vortices_pro,rho_Near_wall_vortices_pro] = ksdensity(Near_wall_vortices_ASL_pro(:,1),...
        linspace(min(Near_wall_vortices_ASL_pro(:,1)), max(Near_wall_vortices_ASL_pro(:,1)), 25),'NumPoints',25);
    [pdfNear_wall_vortices_retro,rho_Near_wall_vortices_retro] = ksdensity(Near_wall_vortices_ASL_retro(:,1),...
        linspace(min(Near_wall_vortices_ASL_retro(:,1)), max(Near_wall_vortices_ASL_retro(:,1)), 25),'NumPoints',25);
    
    pdf_Full_stat_Near_wall_vortices_ASL_pro(z,:) = pdf_rho_Near_wall_vortices_pro;
    rho_Full_stat_Near_wall_vortices_ASL_pro(z,:) = rho_Near_wall_vortices_pro;
    
    pdf_Full_stat_Near_wall_vortices_ASL_retro(z,:) = pdfNear_wall_vortices_retro;
    rho_Full_stat_Near_wall_vortices_ASL_retro(z,:) = rho_Near_wall_vortices_retro;
    z = z+1;
end
%%
m =2;

Near_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2) <= m*VF_PIVWT7.ks,:);
Near_wall_vortices_WT7_pro = Near_wall_vortices_WT7(Near_wall_vortices_WT7(:,3)==-1,:);

Near_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2) <= m*VF_PIVWT10.ks,:);
Near_wall_vortices_WT10_pro = Near_wall_vortices_WT10(Near_wall_vortices_WT10(:,3)==-1,:);

Far_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2) > m*VF_PIVWT7.ks,:);
Far_wall_vortices_WT7_pro = Far_wall_vortices_WT7(Far_wall_vortices_WT7(:,3)==-1,:);

Far_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2) > m*VF_PIVWT10.ks,:);
Far_wall_vortices_WT10_pro = Far_wall_vortices_WT10(Far_wall_vortices_WT10(:,3)==-1,:);

% Define the PDF of the skewed Gaussian distribution
skewGaussPDF_Near = @(x, mu, sigma, alpha) ...
    (1 - x.^2) .* (2 / sigma * normpdf((x - mu) / sigma) .* normcdf(alpha * (x - mu) / sigma));
% skewGaussPDF_Near = @(x, mu, sigma, alpha) ...
%      (2 / sigma * normpdf((x - mu) / sigma) .* normcdf(alpha * (x - mu) / sigma));

% Initial guess for parameters [mu, sigma, alpha]
% initParams = [mean(Near_wall_vortices_WT7_pro(:,1)), std(Near_wall_vortices_WT7_pro(:,1)), 0]; % mu = mean, sigma = std, alpha = skewness
initParams = [mean([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]),...
    std([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]), -0.75]; % mu = mean, sigma = std, alpha = skewness
% Perform curve fitting using 'fit' function (from Curve Fitting Toolbox)
% x_Near = linspace(min(Near_wall_vortices_WT7_pro(:,1)), max(Near_wall_vortices_WT7_pro(:,1)), 100);
x_Near = linspace(min([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]),...
    max([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]), 100);
% y_Near = histcounts(Near_wall_vortices_WT7_pro(:,1), 'Normalization', 'pdf', 'BinEdges',...
%     linspace(min(Near_wall_vortices_WT7_pro(:,1)), max(Near_wall_vortices_WT7_pro(:,1)), 100));
y_Near = histcounts([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)], 'Normalization', 'pdf', 'BinEdges',...
    linspace(min([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]),...
    max([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]), 100));
binCenters_Near = (x_Near(1:end-1) + x_Near(2:end)) / 2; % Bin centers



% Define the PDF of the skewed Gaussian distribution with normalization

% Compute the normalization factor using numerical integration
normalizationFactor = integral(@(x) skewGaussPDF_Near(x, initParams(1), initParams(2), initParams(3)), -1, 1);

% Normalize the function
skewGaussPDF_Near = @(x, mu, sigma, alpha) ...
    (1 - x.^2) .* (2 ./ sigma .* normpdf((x - mu) ./ sigma) .* normcdf(alpha .* (x - mu) ./ sigma)) ./ normalizationFactor;
% skewGaussPDF_Near = @(x, mu, sigma, alpha) ...
%      (2 ./ sigma .* normpdf((x - mu) ./ sigma) .* normcdf(alpha .* (x - mu) ./ sigma)) ./ normalizationFactor;
%Define fit options
fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', initParams, ...
    'Lower', [-Inf, 0, -Inf], 'Upper', [Inf, Inf, Inf]);

fitType = fittype(@(mu, sigma, alpha, x) ...
    skewGaussPDF_Near(x, mu, sigma, alpha), 'independent', 'x', 'options', fitOptions);

% Fit the skewed Gaussian to the histogram data
[fitResult_Near, ~] = fit(binCenters_Near.', y_Near.', fitType);

% Display the fit parameters
disp(fitResult_Near);
skewGaussPDF_Near = @(x) ...
    (1 - x.^2) .* (2 ./ fitResult_Near.sigma .* normpdf((x - fitResult_Near.mu) ./ fitResult_Near.sigma) .* normcdf(fitResult_Near.alpha .* (x - fitResult_Near.mu) ./ fitResult_Near.sigma));
normalizationFactor = integral(@(x) skewGaussPDF_Near(x), -1, 1);
skewGaussPDF_Near = @(x) ...
    (1 - x.^2) .* (2 ./ fitResult_Near.sigma .* normpdf((x - fitResult_Near.mu) ./ fitResult_Near.sigma) .* normcdf(fitResult_Near.alpha .* (x - fitResult_Near.mu) ./ fitResult_Near.sigma))./normalizationFactor;




% Define the PDF of the skewed Gaussian distribution
skewGaussPDF_Far = @(x, mu, sigma, alpha) ...
    (1 - x.^2) .* (2 / sigma * normpdf((x - mu) / sigma) .* normcdf(alpha * (x - mu) / sigma));
% skewGaussPDF_Near = @(x, mu, sigma, alpha) ...
%      (2 ./ sigma .* normpdf((x - mu) ./ sigma) .* normcdf(alpha .* (x - mu) ./ sigma)) ./ normalizationFactor;

% Initial guess for parameters [mu, sigma, alpha]
% initParams = [mean(Far_wall_vortices_WT7_pro(:,1)), std(Far_wall_vortices_WT7_pro(:,1)), 0]; % mu = mean, sigma = std, alpha = skewness
initParams = [mean([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]),...
    std([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]), -100]; % mu = mean, sigma = std, alpha = skewness

% Perform curve fitting using 'fit' function (from Curve Fitting Toolbox)
% x_Far = linspace(min(Far_wall_vortices_WT7_pro(:,1)), max(Far_wall_vortices_WT7_pro(:,1)), 100);
x_Far = linspace(min([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]),...
    max([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]), 100);
% y_Far = histcounts(Far_wall_vortices_WT7_pro(:,1), 'Normalization', 'pdf', 'BinEdges',...
%     linspace(min(Far_wall_vortices_WT7_pro(:,1)), max(Far_wall_vortices_WT7_pro(:,1)), 100));
y_Far = histcounts([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)], 'Normalization', 'pdf', 'BinEdges',...
    linspace(min([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]),...
    max([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]), 100));
binCenters_Far = (x_Far(1:end-1) + x_Far(2:end)) / 2; % Bin centers


normalizationFactor = integral(@(x) skewGaussPDF_Far(x, initParams(1), initParams(2), initParams(3)), -1, 1);

% Normalize the function
skewGaussPDF_Far = @(x, mu, sigma, alpha) ...
    (1 - x.^2) .* (2 ./ sigma .* normpdf((x - mu) ./ sigma) .* normcdf(alpha .* (x - mu) ./ sigma)) ./ normalizationFactor;
% skewGaussPDF_Far = @(x, mu, sigma, alpha) ...
%     (2 ./ sigma .* normpdf((x - mu) ./ sigma) .* normcdf(alpha .* (x - mu) ./ sigma)) ./ normalizationFactor;
%Define fit options
fitOptions = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', initParams, ...
    'Lower', [-Inf, 0, -Inf], 'Upper', [Inf, Inf, Inf]);

%Create a fit type for the skewed Gaussian PDF
fitType = fittype(@(mu, sigma, alpha, x) ...
    skewGaussPDF_Far(x, mu, sigma, alpha), 'independent', 'x', 'options', fitOptions);

% Fit the skewed Gaussian to the histogram data
[fitResult_Far, ~] = fit(binCenters_Far.', y_Far.', fitType);

% Display the fit parameters
disp(fitResult_Far);
skewGaussPDF_Far = @(x) ...
    (1 - x.^2) .* (2 ./ fitResult_Far.sigma .* normpdf((x - fitResult_Far.mu) ./ fitResult_Far.sigma) .* normcdf(fitResult_Far.alpha .* (x - fitResult_Far.mu) ./ fitResult_Far.sigma));
normalizationFactor = integral(@(x) skewGaussPDF_Far(x), -1, 1);
skewGaussPDF_Far = @(x) ...
    (1 - x.^2) .* (2 ./ fitResult_Far.sigma .* normpdf((x - fitResult_Far.mu) ./ fitResult_Far.sigma) .* normcdf(fitResult_Far.alpha .* (x - fitResult_Far.mu) ./ fitResult_Far.sigma))./normalizationFactor;




[pdf_rho_Near_wall_vortices_WT7_WT10_pro,rho_Near_wall_vortices_WT7_WT10_pro] = ksdensity([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)],...
    linspace(min([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]),...
    max([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)]), 25),'NumPoints',25);
[pdf_rho_Near_wall_vortices_WT7_retro,rho_Near_wall_vortices_WT7_retro] = ksdensity(Near_wall_vortices_WT7_retro(:,1),...
        linspace(min(Near_wall_vortices_WT7_retro(:,1)), max(Near_wall_vortices_WT7_retro(:,1)), 25),'NumPoints',25);

    
[pdf_rho_Far_wall_vortices_WT7_WT10_pro,rho_Far_wall_vortices_WT7_WT10_pro] = ksdensity([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)],...
    linspace(min([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]),...
    max([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)]), 25),'NumPoints',25);
[pdf_rho_Far_wall_vortices_WT7_retro,rho_Far_wall_vortices_WT7_retro] = ksdensity(Far_wall_vortices_WT7_retro(:,1),...
        linspace(min(Far_wall_vortices_WT7_retro(:,1)), max(Far_wall_vortices_WT7_retro(:,1)), 25),'NumPoints',25);
figure
subplot(1,2,1)
plot(rho_Near_wall_vortices_WT7_WT10_pro,pdf_rho_Near_wall_vortices_WT7_WT10_pro,...
    'LineStyle','none','color','b','Marker','^');
hold on
plot(rho_Near_wall_vortices_WT7_retro,pdf_rho_Near_wall_vortices_WT7_retro,...
    'LineStyle','none','color','r','Marker','*');
plot(linspace(-1, 1, 100), skewGaussPDF_Near(linspace(-1, 1, 100)),...
    '-k', 'LineWidth', 2);

set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('$\rho$','Interpreter','Latex','FontSize',14);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
legend('Near wall Prograge','Near wall Retrograde','Fitted','Interpreter','latex','FontSize',9)
subplot(1,2,2)
plot(rho_Far_wall_vortices_WT7_WT10_pro,pdf_rho_Far_wall_vortices_WT7_WT10_pro,...
    'LineStyle','none','color','b','Marker','^');
hold on
plot(rho_Far_wall_vortices_WT7_retro,pdf_rho_Far_wall_vortices_WT7_retro,...
    'LineStyle','none','color','r','Marker','*');
plot(linspace(-1, 1, 100), skewGaussPDF_Far(linspace(-1, 1, 100)),...
    '-k', 'LineWidth', 2);
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('$\rho$','Interpreter','Latex','FontSize',14);
ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
legend('Far wall Prograge','Far wall Retrograde','Fitted',...
    'Interpreter','latex','FontSize',9)

%% Hist analysis of pro vs retro r_omega/\lambda for WT7 less than and more than 2ks
m=2;

Near_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*VF_PIVWT7.ks);
Near_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Near_wall_pro_WT7] = min(Near_wall_diffpro, [], 2);

Near_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*VF_PIVWT10.ks);
Near_wall_diffpro = abs(WT10VorX.z_pro(WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Near_wall_pro_WT10] = min(Near_wall_diffpro, [], 2);

Far_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro>m*VF_PIVWT7.ks);
Far_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro>m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Far_wall_pro_WT7] = min(Far_wall_diffpro, [], 2);

Far_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro > m*VF_PIVWT10.ks);
Far_wall_diffpro = abs(WT10VorX.z_pro(WT10VorX.z_pro > m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Far_wall_pro_WT10] = min(Far_wall_diffpro, [], 2);

data_WT7_near_pro = double(0.5*Near_wall_d_pro_WT7./WT7VorX_scales.lambda(result_Near_wall_pro_WT7));
data_WT10_near_pro = double(0.5*Near_wall_d_pro_WT10./WT10VorX_scales.lambda(result_Near_wall_pro_WT10));
data_WT7_Far_pro = double(0.5*Far_wall_d_pro_WT7./WT7VorX_scales.lambda(result_Far_wall_pro_WT7));
data_WT10_Far_pro = double(0.5*Far_wall_d_pro_WT10./WT10VorX_scales.lambda(result_Far_wall_pro_WT10));

Near_wall_d_ret_WT7 = WT7VorX.d_ret (WT7VorX.z_ret>=0 & WT7VorX.z_ret <= m*VF_PIVWT7.ks);
Near_wall_diffret = abs(WT7VorX.z_ret(WT7VorX.z_ret>=0 & WT7VorX.z_ret <= m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Near_wall_ret_WT7] = min(Near_wall_diffret, [], 2);

Near_wall_d_ret_WT10 = WT10VorX.d_ret (WT10VorX.z_ret>=0 & WT10VorX.z_ret <= m*VF_PIVWT10.ks);
Near_wall_diffret = abs(WT10VorX.z_ret(WT10VorX.z_ret>=0 & WT10VorX.z_ret<= m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Near_wall_ret_WT10] = min(Near_wall_diffret, [], 2);

Far_wall_d_ret_WT7 = WT7VorX.d_ret (WT7VorX.z_ret>m*VF_PIVWT7.ks);
Far_wall_diffret = abs(WT7VorX.z_ret(WT7VorX.z_ret>m*VF_PIVWT7.ks)...
    - WT7VorX_scales.z.');
[~, result_Far_wall_ret_WT7] = min(Far_wall_diffret, [], 2);

Far_wall_d_ret_WT10 = WT10VorX.d_ret (WT10VorX.z_ret>m*VF_PIVWT10.ks);
Far_wall_diffret = abs(WT10VorX.z_ret(WT10VorX.z_ret>m*VF_PIVWT10.ks)...
    - WT10VorX_scales.z.');
[~, result_Far_wall_ret_WT10] = min(Far_wall_diffret, [], 2);


data_WT7_near_ret = double(0.5*Near_wall_d_ret_WT7./WT7VorX_scales.lambda(result_Near_wall_ret_WT7));
data_WT10_near_ret = double(0.5*Near_wall_d_ret_WT10./WT10VorX_scales.lambda(result_Near_wall_ret_WT10));
data_WT7_Far_ret = double(0.5*Far_wall_d_ret_WT7./WT7VorX_scales.lambda(result_Far_wall_ret_WT7));
data_WT10_Far_ret = double(0.5*Far_wall_d_ret_WT10./WT10VorX_scales.lambda(result_Far_wall_ret_WT10));

Y_edges = logspace(log10(3e-2), log10(3e0), 25);
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end));

figure
set(gcf,'Position',[689,419,806,545])
axes('Position',[0.413151364764265,0.118790496760259,0.241818730910886,0.741759961955337])
[N_near_pro,~] = histcounts([data_WT7_near_pro;data_WT10_near_pro],Y_edges,'Normalization','pdf');
plot(Y_centers,N_near_pro,'r*','MarkerSize',9,'LineWidth',1.5)
xlim([3*10^-2 3*10^0])
ylim([1e-5 2e1])
hold on
[N_near_ret,~] = histcounts([data_WT7_near_ret;data_WT10_near_ret],Y_edges,'Normalization','pdf');
plot(Y_centers,N_near_ret,'bO','MarkerSize',9,'LineWidth',1.5)

plot(x_values_near,fitted_pdf_near,'Color',[0.00,0.00,0.00],'LineWidth',3,...
    'LineStyle','-')

[N_Far_pro,~] = histcounts([data_WT7_Far_pro;data_WT10_Far_pro],Y_edges,'Normalization','pdf');
plot(Y_centers,N_Far_pro,'r^','MarkerSize',9,'LineWidth',1.5)


[N_Far_ret,~] = histcounts([data_WT7_Far_ret;data_WT10_Far_ret],Y_edges,'Normalization','pdf');
plot(Y_centers,N_Far_ret,'bs','MarkerSize',11,'LineWidth',1.5)


plot(x_values_Far,fitted_pdf_Far,'Color',[0.00,1.00,0.00],'LineWidth',3,...
    'LineStyle','-.')

set(gca,'XScale','log','YScale','log','XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on',...
    'TickLabelInterpreter','latex','FontSize',15)
legend('$\mathrm{z}/\mathrm{k}_{\mathrm{s}}\le2$ Prograde','$\mathrm{z}/\mathrm{k}_{\mathrm{s}}\le2$ Retrograde','Fitted  $\mathrm{z}/\mathrm{k}_{s}<2$',...
    '$\mathrm{z}/\mathrm{k}_{\mathrm{s}}\ge2$ Prograde','$\mathrm{z}/\mathrm{k}_{\mathrm{s}}\ge2$ Retrograde','Fitted  $\mathrm{z}/\mathrm{k}_{s}\geq2$',...
    'Interpreter', 'latex','NumColumns',2,'position',[0.006203473945409,0.869241734792314,0.988833746898263,0.128073398301361])
xlabel('$\mathrm{r}_{\omega}/\lambda_{T}$','Interpreter','Latex','FontSize',17)
ylabel('p.d.f.','Interpreter','Latex','FontSize',17)
annotation(gcf,'textbox',...
    [0.602736972704715,0.787155963302753,0.046146401985112,0.06605504587156],...
    'String',{'(b)'},...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Hist analysis of pro vs retro u_omega/u_tau for WT7 less than and more than 2ks
m=2;

Near_wall_U_omega_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*VF_PIVWT7.ks)/VF_PIVWT7.u_tau;
Near_wall_U_omega_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*VF_PIVWT10.ks)/VF_PIVWT10.u_tau;
Far_wall_U_omega_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro > m*VF_PIVWT7.ks)/VF_PIVWT7.u_tau;
Far_wall_U_omega_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro > m*VF_PIVWT10.ks)/VF_PIVWT10.u_tau;

data_WT7_near_pro = double(-Near_wall_U_omega_pro_WT7);
data_WT10_near_pro = double(-Near_wall_U_omega_pro_WT10);
data_WT7_Far_pro = double(-Far_wall_U_omega_pro_WT7);
data_WT10_Far_pro = double(-Far_wall_U_omega_pro_WT10);

Near_wall_U_omega_ret_WT7 = WT7VorX.dU_ret (WT7VorX.z_ret>=0 & WT7VorX.z_ret <= m*VF_PIVWT7.ks)/VF_PIVWT7.u_tau;
Near_wall_U_omega_ret_WT10 = WT10VorX.dU_ret (WT10VorX.z_ret>=0 & WT10VorX.z_ret <= m*VF_PIVWT10.ks)/VF_PIVWT10.u_tau;
Far_wall_U_omega_ret_WT7 = WT7VorX.dU_ret (WT7VorX.z_ret > m*VF_PIVWT7.ks)/VF_PIVWT7.u_tau;
Far_wall_U_omega_ret_WT10 = WT10VorX.dU_ret (WT10VorX.z_ret > m*VF_PIVWT10.ks)/VF_PIVWT10.u_tau;


data_WT7_near_ret = double(Near_wall_U_omega_ret_WT7);
data_WT10_near_ret = double(Near_wall_U_omega_ret_WT10);
data_WT7_Far_ret = double(Far_wall_U_omega_ret_WT7);
data_WT10_Far_ret = double(Far_wall_U_omega_ret_WT10);



axes('Position',[0.750620347394541,0.118790496760259,0.241818730910886,0.741759961955337])
[N_near_pro,edges] = histcounts([data_WT7_near_pro;data_WT10_near_pro],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_near_pro,'r*','MarkerSize',9,'LineWidth',1.5)
xlim([-0.5 6])
ylim([1e-6 1e1])
hold on
[N_near_ret,edges] = histcounts([data_WT7_near_ret;data_WT10_near_ret],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_near_ret,'bO','MarkerSize',9,'LineWidth',1.5)


% semilogy((0.5:0.1:4.5),5.43*exp(-(0.5:0.1:4.5)./(1*0.5)),'Color',[0.00,0.00,0.00],'LineWidth',3)
semilogy(x_values_near_u_omega,fitted_pdf_near_u_omega,'Color',[0.00,0.00,0.00],'LineWidth',3)

[N_Far_pro,edges] = histcounts([data_WT7_Far_pro;data_WT10_Far_pro],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_Far_pro,'r^','MarkerSize',9,'LineWidth',1.5)


[N_Far_ret,edges] = histcounts([data_WT7_Far_ret;data_WT10_Far_ret],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_Far_ret,'bs','MarkerSize',11,'LineWidth',1.5)

semilogy(x_values_Far_u_omega,fitted_pdf_Far_u_omega,'Color',[0.00,1.00,0.00],'LineWidth',3,...
    'LineStyle','-.')
% semilogy((0.5:0.1:4.5),10.43*exp(-(0.5:0.1:4.5)./(1*0.37)),'Color',[0.00,1.00,0.00],'LineWidth',3,...
%     'LineStyle','-.')

set(gca,'XScale','linear','YScale','log','XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on',...
    'TickLabelInterpreter','latex','FontSize',15)
xlabel('$\mathrm{u}_{\omega}/u_{\tau}$','Interpreter','Latex','FontSize',17)
ylabel('p.d.f.','Interpreter','Latex','FontSize',17)
annotation(gcf,'textbox',...
    [0.940205955334988,0.787155963302753,0.046146401985112,0.06605504587156],...
    'String',{'(c)'},...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Hist analysis of pro vs retro \rho for WT7 less than and more than 2ks

Near_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2) <= m*VF_PIVWT7.ks,:);
Near_wall_vortices_WT7_pro = Near_wall_vortices_WT7(Near_wall_vortices_WT7(:,3)==-1,:);
Near_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2) <= m*VF_PIVWT10.ks,:);
Near_wall_vortices_WT10_pro = Near_wall_vortices_WT10(Near_wall_vortices_WT10(:,3)==-1,:);

Near_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2) <= m*VF_PIVWT7.ks,:);
Near_wall_vortices_WT7_ret = Near_wall_vortices_WT7(Near_wall_vortices_WT7(:,3)==1,:);
Near_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2) <= m*VF_PIVWT10.ks,:);
Near_wall_vortices_WT10_ret = Near_wall_vortices_WT10(Near_wall_vortices_WT10(:,3)==1,:);

Far_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2) > m*VF_PIVWT7.ks,:);
Far_wall_vortices_WT7_pro = Far_wall_vortices_WT7(Far_wall_vortices_WT7(:,3)==-1,:);
Far_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2) > m*VF_PIVWT10.ks,:);
Far_wall_vortices_WT10_pro = Far_wall_vortices_WT10(Far_wall_vortices_WT10(:,3)==-1,:);

Far_wall_vortices_WT7 = corr_coef_local_WT7 (corr_coef_local_WT7(:,2) > m*VF_PIVWT7.ks,:);
Far_wall_vortices_WT7_ret = Far_wall_vortices_WT7(Far_wall_vortices_WT7(:,3)==1,:);
Far_wall_vortices_WT10 = corr_coef_local_WT10 (corr_coef_local_WT10(:,2) > m*VF_PIVWT10.ks,:);
Far_wall_vortices_WT10_ret = Far_wall_vortices_WT10(Far_wall_vortices_WT10(:,3)==1,:);

axes('Position',[0.076923076923073,0.118790496760259,0.241818730910885,0.741759961955337])
[N_near_pro,edges] = histcounts([Near_wall_vortices_WT7_pro(:,1);Near_wall_vortices_WT10_pro(:,1)],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_near_pro,'r*','MarkerSize',9,'LineWidth',1.5)
xlim([-1 1])
ylim([0 1.3])
hold on
[N_near_ret,edges] = histcounts([Near_wall_vortices_WT7_ret(:,1);Near_wall_vortices_WT10_ret(:,1)],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_near_ret,'bO','MarkerSize',9,'LineWidth',1.5)

plot(linspace(-1, 1, 100), skewGaussPDF_Near(linspace(-1, 1, 100)),'Color',[0.00,0.00,0.00],'LineWidth',3,...
    'LineStyle','-')

[N_Far_pro,edges] = histcounts([Far_wall_vortices_WT7_pro(:,1);Far_wall_vortices_WT10_pro(:,1)],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_Far_pro,'r^','MarkerSize',9,'LineWidth',1.5)

[N_Far_ret,edges] = histcounts([Far_wall_vortices_WT7_ret(:,1);Far_wall_vortices_WT10_ret(:,1)],25,'Normalization','pdf');
plot(0.5*(edges(1,1:end-1)+edges(1,2:end)),N_Far_ret,'bs','MarkerSize',11,'LineWidth',1.5)

plot(linspace(-1, 1, 100), skewGaussPDF_Far(linspace(-1, 1, 100)),'Color',[0.00,1.00,0.00],'LineWidth',3,...
    'LineStyle','-.')

set(gca,'XScale','linear','YScale','linear','XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on',...
    'TickLabelInterpreter','latex','FontSize',15)
xlabel('$\rho_{\mathrm{u,w}}^{\mathrm{VorX}}$','Interpreter','Latex','FontSize',17)
ylabel('p.d.f.','Interpreter','Latex','FontSize',17)
annotation(gcf,'textbox',...
    [0.266508684863524,0.787155963302753,0.046146401985112,0.06605504587156],...
    'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Generating random numbers Near wall

% Step 1: Define the range and resolution for numerical integration
x_vals = linspace(-1, 1, 1000);  % High resolution for accuracy

% Step 2: Calculate the PDF over the range
pdf_vals = skewGaussPDF_Near(x_vals);

% Step 3: Compute the CDF using cumulative trapezoidal integration
cdf_vals = cumtrapz(x_vals, pdf_vals);
cdf_vals = cdf_vals ./ max(cdf_vals);  % Normalize to ensure CDF ends at 1

% Step 4: Invert the CDF using interpolation
inverseCDF = @(u) interp1(cdf_vals, x_vals, u, 'linear', 'extrap');

% Step 5: Generate uniform random numbers
nSamples = 1e7;  % Number of samples
u = rand(nSamples, 1);  % Uniform random numbers in [0, 1]

% Step 6: Generate samples using the inverse CDF
samples_Near = inverseCDF(u);

figure;
histogram(samples_Near, 100, 'Normalization', 'pdf'); % Histogram of generated samples
hold on;

% Generate x values for plotting the PDF
x = linspace(min(samples_Near), max(samples_Near), 100);
plot(x, skewGaussPDF_Near(x), '-r', 'LineWidth', 2); % Fitted PDF

% Add plot labels and legend
title('Generated Samples from Fitted Skewed Gaussian');
xlabel('Data Values');
ylabel('Probability Density');
legend('Generated Data', 'Fitted PDF');
hold off;

%% Generating random numbers Far wall

% Step 1: Define the range and resolution for numerical integration
x_vals = linspace(-1, 1, 1000);  % High resolution for accuracy

% Step 2: Calculate the PDF over the range
pdf_vals = skewGaussPDF_Far(x_vals);

% Step 3: Compute the CDF using cumulative trapezoidal integration
cdf_vals = cumtrapz(x_vals, pdf_vals);
cdf_vals = cdf_vals ./ max(cdf_vals);  % Normalize to ensure CDF ends at 1

% Step 4: Invert the CDF using interpolation
inverseCDF = @(u) interp1(cdf_vals, x_vals, u, 'linear', 'extrap');

% Step 5: Generate uniform random numbers
nSamples = 1e7;  % Number of samples
u = rand(nSamples, 1);  % Uniform random numbers in [0, 1]

% Step 6: Generate samples using the inverse CDF
samples_Far = inverseCDF(u);

figure;
histogram(samples_Far, 100, 'Normalization', 'pdf'); % Histogram of generated samples
hold on;

% Generate x values for plotting the PDF
x = linspace(min(samples_Far), max(samples_Far), 100);
plot(x, skewGaussPDF_Far(x), '-r', 'LineWidth', 2); % Fitted PDF

% Add plot labels and legend
title('Generated Samples from Fitted Skewed Gaussian');
xlabel('Data Values');
ylabel('Probability Density');
legend('Generated Data', 'Fitted PDF');
hold off;
%%

figure;
histogram(samples_Near, 'Normalization', 'pdf'); 
hold on;
[f, xi] = ksdensity(samples_Near); 
plot(xi, f, 'r-', 'LineWidth', 2); 


figure;
histogram(samples_Far, 'Normalization', 'pdf'); 
hold on;
[f, xi] = ksdensity(samples_Far); 
plot(xi, f, 'r-', 'LineWidth', 2); 
%%
save('samples_Near.mat','samples_Near')
save('samples_Far.mat','samples_Far')

load('samples_Near.mat')
load('samples_Far.mat')
