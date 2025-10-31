Data1 = load('VF_GenWT7(Beforepassingkernelres32makima).mat');
VF_GenWT7 = Data1.VF_GenWT7;
clear('Data1')

% Data2 = load('VF_PIVWT7.mat');
% VF_PIVWT7 = Data2.VF_PIVWT7;
% clear('Data2')

Data3 = load('VF_HotWT7.mat');
VF_HotWT7 = Data3.VF_HotWT7;
clear('Data3')

Data4 = load('VF_GenWT10(Beforepassingkernelres25makima).mat');
VF_GenWT10 = Data4.VF_GenWT10;
clear('Data4')

% Data5 = load('VF_PIVWT10.mat');
% VF_PIVWT10 = Data5.VF_PIVWT10;
% clear('Data5')

Data6 = load('VF_HotWT10.mat');
VF_HotWT10 = Data6.VF_HotWT10;
clear('Data6')

Data7 = load('VF_GenASL(Beforepassingkernelres20makima).mat');
VF_GenASL = Data7.VF_GenASL;
clear('Data7')

Data8 = load('VF_SLPIVASL.mat');
VF_SLPIVASL = Data8.VF_SLPIVASL;
clear('Data8')

Data9 = load('VF_SonicASL.mat');
VF_SonicASL = Data9.VF_SonicASL;
clear('Data9')

%% Calculating lambda_T
lambda_T_WT7 = zeros(size(VF_HotWT7.uprime, 1), 1);
eta_WT7 = zeros(size(VF_HotWT7.uprime, 1), 1);

for i = 1:size(VF_HotWT7.uprime, 1)

    current_row_epsilon = mean([VF_HotWT7.epsilon_str{i, :}{:}],2);
    
    lambda_T_WT7(i) = rms(VF_HotWT7.uprime(i,:),2)*sqrt(15*VF_HotWT7.nu/current_row_epsilon);
    eta_WT7(i)= (VF_HotWT7.nu^3/current_row_epsilon)^0.25;
end


lambda_T_WT10 = zeros(size(VF_HotWT10.uprime, 1), 1);
eta_WT10 = zeros(size(VF_HotWT10.uprime, 1), 1);

for i = 1:size(VF_HotWT10.uprime, 1)

    current_row_epsilon = mean([VF_HotWT10.epsilon_str{i, :}{:}],2);
    
    lambda_T_WT10(i) = rms(VF_HotWT10.uprime(i,:),2)*sqrt(15*VF_HotWT10.nu/current_row_epsilon);
    eta_WT10(i)= (VF_HotWT10.nu^3/current_row_epsilon)^0.25;
end


current_row_epsilon_SonicASL = mean([VF_SonicASL.epsilon_str{1, :}{:}],2);
lambda_T_ASL = rms(VF_SonicASL.uprime(:,1),1)*sqrt(15*VF_SonicASL.nu/current_row_epsilon_SonicASL);
eta_ASL = (VF_SonicASL.nu^3/current_row_epsilon_SonicASL)^0.25;
% lambda_T_ASL = rms(VF_SLPIVASL.uprime(7,25,:),3)*sqrt(15*VF_SonicASL.nu/current_row_epsilon_SonicASL)

%% Choping the signal into lambda_T

T_domain_WT7 = (0:1:size(VF_HotWT7.u,2)-1)*1/VF_HotWT7.fs;
X_domain_WT7 = mean(VF_HotWT7.u,2)*T_domain_WT7;

T_domain_WT10 = (0:1:size(VF_HotWT10.u,2)-1)*1/VF_HotWT10.fs;
X_domain_WT10 = mean(VF_HotWT10.u,2)*T_domain_WT10;

T_domain_ASL = (0:1:size(VF_SLPIVASL.u,3)-1)*1/VF_SLPIVASL.fs;
T_domain_ASL = reshape(T_domain_ASL,1,1,[]);
X_domain_ASL = mean(VF_SLPIVASL.u,3).*T_domain_ASL;

VF_HotWT7uprimeshort = zeros(size(VF_HotWT7.u));
VF_HotWT7wprimeshort = zeros(size(VF_HotWT7.w));
VF_HotWT10uprimeshort = zeros(size(VF_HotWT10.u));
VF_HotWT10wprimeshort = zeros(size(VF_HotWT10.w));
VF_SLPIVASLuprimeshort = zeros(size(VF_SLPIVASL.u));
VF_SLPIVASLwprimeshort = zeros(size(VF_SLPIVASL.w));
C = zeros();
% figure
% set(gcf,'Position',[814,534,806,379])
% axes('Position',[0.064516129032258,0.12664907651715,0.915632754342432,0.798350923482849])

l0WT7 = 0.7;%med = 0.3, long = 0.7, long-long = 1.2
l0WT10 = 0.7;%med = 0.3,long = 0.7, long-long = 1.2
l0ASL = 0.72;


c0WT7 = 3;
c0WT10 = 3;
c0ASL = 3;

for z = 1: size(VF_HotWT7.z,1)
    
%     ratio_WT7= X_domain_WT7(z,:) / (c0WT7*lambda_T_WT7(z,1));
    ratio_WT7= X_domain_WT7(z,:) / l0WT7;
    unitindex_WT7 = find(ratio_WT7 >= 1, 1, 'first');
    unitindex_WT7 = unitindex_WT7 - 1;
    
    for s = 1 : floor(size(VF_HotWT7.u(z,:),2)/unitindex_WT7)
        VF_HotWT7uprimeshort(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2)+1e-6*rand(1,unitindex_WT7);
        
        VF_HotWT7wprimeshort(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2)+1e-6*rand(1,unitindex_WT7);
        
        
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/0.3,...
%             VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
%         hold on
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/0.3,...
%             VF_HotWT7uprimeshort(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
%         
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/0.3,...
%             VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-mean(VF_HotWT7.u(z,:),2),'linewidth',2)
%         legend(sprintf('u[m/s],H-W(m1) z$/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
%             sprintf('$\\mathrm{u}^{\\prime}_{\\mathrm{Local}}\\mathrm{[m/s]}$,H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
%             sprintf('$\\mathrm{u}^{\\prime}_{\\mathrm{Global}}\\mathrm{[m/s]}$,H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
%             'Interpreter','latex','FontSize',10,'Position',...
%             [0.69418761794644,0.810046174733969,0.271111721158469,0.096306066085292],...
%             'Numcolumns',1,'Orientation','vertical','color','none');
%         xlabel('x/L','Interpreter','latex')
%         ylabel('$\mathrm{u},\mathrm{u}^{\prime}$[m/s]','Interpreter','latex')
%         set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on','YGrid','on')
        
        R = corrcoef(VF_HotWT7uprimeshort(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)',...
            VF_HotWT7wprimeshort(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)');
        C(z,s) = R(1,2);
    end
    %Add ks density
    
end

cols_with_zeros_uprime_WT7 = any(VF_HotWT7uprimeshort == 0, 1);
VF_HotWT7uprimeshort = VF_HotWT7uprimeshort(:, ~cols_with_zeros_uprime_WT7);

cols_with_zeros_wprime_WT7 = any(VF_HotWT7wprimeshort == 0, 1);
VF_HotWT7wprimeshort = VF_HotWT7wprimeshort(:, ~cols_with_zeros_wprime_WT7);

for z=1:size(VF_HotWT7.z,1)
    VF_HotWT7.z(z)/VF_HotWT7.delta
    R = corrcoef(VF_HotWT7uprimeshort(z,:),VF_HotWT7wprimeshort(z,:))
    
end

for z = 1: size(VF_HotWT10.z,1)
    
%     ratio_WT10 = X_domain_WT10(z,:) / (c0WT10*lambda_T_WT10(z,1));
    ratio_WT10 = X_domain_WT10(z,:) / l0WT10;
    unitindex_WT10 = find(ratio_WT10 >= 1, 1, 'first');
    unitindex_WT10 = unitindex_WT10 - 1;
    
    for s = 1 : floor(size(VF_HotWT10.u(z,:),2)/unitindex_WT10)
        VF_HotWT10uprimeshort(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2)+1e-6*rand(1,unitindex_WT10);
        
        VF_HotWT10wprimeshort(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2)+1e-6*rand(1,unitindex_WT10);        
    end
    
end

cols_with_zeros_uprime_WT10 = any(VF_HotWT10uprimeshort == 0, 1);
VF_HotWT10uprimeshort = VF_HotWT10uprimeshort(:, ~cols_with_zeros_uprime_WT10);

cols_with_zeros_wprime_WT10 = any(VF_HotWT10wprimeshort == 0, 1);
VF_HotWT10wprimeshort = VF_HotWT10wprimeshort(:, ~cols_with_zeros_wprime_WT10);

for z=1:size(VF_HotWT10.z,1)
    
    VF_HotWT10.z(z)/VF_HotWT10.delta
    R = corrcoef(VF_HotWT10uprimeshort(z,:),VF_HotWT10wprimeshort(z,:))
    
end

for z = 1: size(VF_SLPIVASL.z,1)
    
%     ratio_ASL = X_domain_ASL(z,:,:) ./ lambda_T_ASL(1,1);
    ratio_ASL = X_domain_ASL(z,:,:) ./ l0ASL;
    ratio_ASL_squeezed = squeeze(ratio_ASL); % Convert to 52 × 80000 matrix

    % Initialize the result vector with NaN (in case no values > 1 are found)
    unitindex_ASL = NaN(1, size(ratio_ASL_squeezed,1));

    % Loop through each column to find the first index where value > 1
    for col = 1:size(ratio_ASL_squeezed,1)
        idx = find(ratio_ASL_squeezed(col, :) >= 1, 1, 'first'); % Find first occurrence
        if ~isempty(idx)
            unitindex_ASL(col) = idx -1 ; % Store the index if found
        end
    end


    for j = 1:size(unitindex_ASL,2)
        for s = 1 : floor(size(VF_SLPIVASL.u(z,j,:),3)/unitindex_ASL(1,j))
            VF_SLPIVASLuprimeshort(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)) = ...
                VF_SLPIVASL.u(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j))-...
                mean(VF_SLPIVASL.u(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)),3)+1e-6*rand(1,1,unitindex_ASL(1,j));

            VF_SLPIVASLwprimeshort(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)) = ...
                VF_SLPIVASL.w(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j))-...
                mean(VF_SLPIVASL.w(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)),3)+1e-6*rand(1,1,unitindex_ASL(1,j));
            
        end
    end
end

cols_with_zeros_uprime_ASL = squeeze(any(VF_SLPIVASLuprimeshort == 0, [1 2]));
VF_SLPIVASLuprimeshort = VF_SLPIVASLuprimeshort(:,:, ~cols_with_zeros_uprime_ASL);

cols_with_zeros_wprime_ASL = squeeze(any(VF_SLPIVASLwprimeshort == 0, [1 2]));
VF_SLPIVASLwprimeshort = VF_SLPIVASLwprimeshort(:,:, ~cols_with_zeros_wprime_ASL);


%% Analyzing the distriburtion

mean_uprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
std_uprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
mean_wprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
std_wprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);


for z = 1: size(VF_HotWT7.z,1)
    
    
    
    [pdf_uprime_shorttime_WT7,uprime_shorttime_WT7] = ksdensity(VF_HotWT7uprimeshort(z,:)/VF_HotWT7.u_tau,...
    linspace(-6, 6, 100));

    [pdf_uprime_globaltime_WT7,uprime_globaltime_WT7] = ksdensity((VF_HotWT7.u(z,:)-mean(VF_HotWT7.u(z,:),2))/VF_HotWT7.u_tau,...
    linspace(-6, 6, 100));

    [pdf_wprime_shorttime_WT7,wprime_shorttime_WT7] = ksdensity(VF_HotWT7wprimeshort(z,:),...
    linspace(min(VF_HotWT7wprimeshort(z,:)), max(VF_HotWT7wprimeshort(z,:)), 50),'NumPoints',50);


    mean_uprime_shorttime_HotWT7(z,1) = mean(VF_HotWT7uprimeshort(z,:),2);
    std_uprime_shorttime_HotWT7(z,1) = std(VF_HotWT7uprimeshort(z,:),0,2);
    
    mean_wprime_shorttime_HotWT7(z,1) = mean(VF_HotWT7wprimeshort(z,:),2);
    std_wprime_shorttime_HotWT7(z,1) = std(VF_HotWT7wprimeshort(z,:),0,2);
    
    pdf_u_fit = normpdf(linspace(min(VF_HotWT7uprimeshort(z,:)), max(VF_HotWT7uprimeshort(z,:)), 1000), 0, std_uprime_shorttime_HotWT7(z,1));
    pdf_w_fit = normpdf(linspace(min(VF_HotWT7wprimeshort(z,:)), max(VF_HotWT7wprimeshort(z,:)), 1000), 0, std_wprime_shorttime_HotWT7(z,1));
    
%     figure
%     set(gcf,'Position',[622,508,806,394])
%     axes('Position',[0.08560794044665,0.134517766497462,0.40818858560794,0.83502538071066])
%     plot(uprime_shorttime_WT7,pdf_uprime_shorttime_WT7,...
%         'LineStyle','-','color','r','Marker','none','Linewidth',2);
%     hold on
%     plot(uprime_globaltime_WT7,pdf_uprime_globaltime_WT7,...
%         'LineStyle','-.','color','b','Marker','none','Linewidth',2);
% %     plot(linspace(min(VF_HotWT7uprimeshort(z,:)), max(VF_HotWT7uprimeshort(z,:)), 1000), pdf_u_fit, 'r-', 'LineWidth', 1.5); % Normal distribution
%     
%     legend(sprintf('Sub--window mean, H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
%         sprintf('Global mean, H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
%         'Interpreter','latex','FontSize',11,'Position',...
%         [0.094243393010764,0.853855920868091,0.389991769448501,0.100761418173156],...
%         'Numcolumns',1,'Orientation','vertical','color',[1,1,1]);
%     xlabel('$\mathrm{u}^{\prime}/u_{\tau}$','Interpreter','latex')
%     ylabel('p.d.f.','Interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex','FontSize',14,'XGrid','on','YGrid','on')
%     xlim([-5.5 5.5])
%     ylim([0 0.3])
%     axis square
    
%     figure
%     plot(wprime_shorttime_WT7,pdf_wprime_shorttime_WT7,...
%     'LineStyle','none','color','k','Marker','^');
%     hold on
%     plot(linspace(min(VF_HotWT7wprimeshort(z,:)), max(VF_HotWT7wprimeshort(z,:)), 1000), pdf_w_fit, 'k-', 'LineWidth', 1.5); % Normal distribution
%     
end


mean_uprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_uprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
mean_wprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_wprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);

for z = 1: size(VF_HotWT10.z,1)
    
    
    
    [pdf_uprime_shorttime_WT10,uprime_shorttime_WT10] = ksdensity(VF_HotWT10uprimeshort(z,:)/VF_HotWT10.u_tau,...
    linspace(-6, 6, 100));

    [pdf_uprime_globaltime_WT10,uprime_globaltime_WT10] = ksdensity((VF_HotWT10.u(z,:)-mean(VF_HotWT10.u(z,:),2))/VF_HotWT10.u_tau,...
    linspace(-6, 6, 100));

    [pdf_wprime_shorttime_WT10,wprime_shorttime_WT10] = ksdensity(VF_HotWT10wprimeshort(z,:),...
    linspace(min(VF_HotWT10wprimeshort(z,:)), max(VF_HotWT10wprimeshort(z,:)), 50),'NumPoints',50);

    mean_uprime_shorttime_HotWT10(z,1) = mean(VF_HotWT10uprimeshort(z,:),2);
    std_uprime_shorttime_HotWT10(z,1) = std(VF_HotWT10uprimeshort(z,:),0,2);
    
    mean_wprime_shorttime_HotWT10(z,1) = mean(VF_HotWT10wprimeshort(z,:),2);
    std_wprime_shorttime_HotWT10(z,1) = std(VF_HotWT10wprimeshort(z,:),0,2);
    

    pdf_u_fit = normpdf(linspace(min(VF_HotWT10uprimeshort(z,:)), max(VF_HotWT10uprimeshort(z,:)), 1000), 0, std_uprime_shorttime_HotWT10(z,1));
    pdf_w_fit = normpdf(linspace(min(VF_HotWT10wprimeshort(z,:)), max(VF_HotWT10wprimeshort(z,:)), 1000), 0, std_wprime_shorttime_HotWT10(z,1));
    
%     figure
%     set(gcf,'Position',[622,508,806,394])
%     axes('Position',[0.08560794044665,0.134517766497462,0.40818858560794,0.83502538071066])
%     plot(uprime_shorttime_WT10,pdf_uprime_shorttime_WT10,...
%         'LineStyle','-','color','r','Marker','none','Linewidth',2);
%     hold on
%     plot(uprime_globaltime_WT10,pdf_uprime_globaltime_WT10,...
%         'LineStyle','-.','color','b','Marker','none','Linewidth',2);
% %     plot(linspace(min(VF_HotWT10uprimeshort(z,:)), max(VF_HotWT10uprimeshort(z,:)), 1000), pdf_u_fit, 'r-', 'LineWidth', 1.5); % Normal distribution
%     
%     legend(sprintf('Sub--window mean, H-W(m2) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT10.z(z,1)/VF_HotWT10.delta),...
%         sprintf('Global mean, H-W(m2) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT10.z(z,1)/VF_HotWT10.delta),...
%         'Interpreter','latex','FontSize',11,'Position',...
%         [0.094243393010764,0.853855920868091,0.389991769448501,0.100761418173156],...
%         'Numcolumns',1,'Orientation','vertical','color',[1,1,1]);
%     xlabel('$\mathrm{u}^{\prime}/u_{\tau}$','Interpreter','latex')
%     ylabel('p.d.f.','Interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex','FontSize',14,'XGrid','on','YGrid','on')
%     xlim([-5.5 5.5])
%     ylim([0 0.3])
%     axis square


end

mean_uprime_shorttime_SLPIVASL = zeros(size(VF_SLPIVASL.z,1),1);
std_uprime_shorttime_SLPIVASL = zeros(size(VF_SLPIVASL.z,1),1);
mean_wprime_shorttime_SLPIVASL = zeros(size(VF_SLPIVASL.z,1),1);
std_wprime_shorttime_SLPIVASL = zeros(size(VF_SLPIVASL.z,1),1);

for z = 1: size(VF_SLPIVASL.z,1)
    
   [pdf_uprime_shorttime_ASL,uprime_shorttime_ASL] = ksdensity(reshape(VF_SLPIVASLuprimeshort(z,:,:),1,[]),...
    linspace(-7, 7, 100));
    
    [pdf_uprime_globaltime_ASL,uprime_globaltime_ASL] = ksdensity(reshape((VF_SLPIVASL.u(z,:,:)-mean(VF_SLPIVASL.u(z,:,:),3))/VF_SLPIVASL.u_tau,1,[]),...
    linspace(-7, 7, 100));

    [pdf_wprime_shorttime_ASL,wprime_shorttime_ASL] = ksdensity(reshape(VF_SLPIVASLwprimeshort(z,:,:),1,[]),...
    linspace(min(VF_SLPIVASLwprimeshort(z,:),[],'all'), max(VF_SLPIVASLwprimeshort(z,:),[],'all'), 50),'NumPoints',50);

    mean_uprime_shorttime_SLPIVASL(z,1) = mean(mean(VF_SLPIVASLuprimeshort(z,:,:),3),2);
    std_uprime_shorttime_SLPIVASL(z,1) = mean(std(VF_SLPIVASLuprimeshort(z,:,:),0,3),2);
    
    mean_wprime_shorttime_SLPIVASL(z,1) = mean(mean(VF_SLPIVASLwprimeshort(z,:,:),3),2);
    std_wprime_shorttime_SLPIVASL(z,1) = mean(std(VF_SLPIVASLwprimeshort(z,:,:),0,3),2);
    
    pdf_u_fit = normpdf(linspace(min(VF_SLPIVASLuprimeshort(z,:,:),[],'all'),...
        max(VF_SLPIVASLuprimeshort(z,:,:),[],'all'), 1000), 0, std_uprime_shorttime_SLPIVASL(z,1));
    pdf_w_fit = normpdf(linspace(min(VF_SLPIVASLwprimeshort(z,:,:),[],'all'),...
        max(VF_SLPIVASLwprimeshort(z,:,:),[],'all'), 1000), 0, std_wprime_shorttime_SLPIVASL(z,1));
    
%     figure
%     plot(uprime_shorttime_ASL,pdf_uprime_shorttime_ASL,...
%         'LineStyle','-','color','r','Marker','none','Linewidth',2);
%     hold on
%     plot(uprime_globaltime_ASL,pdf_uprime_globaltime_ASL,...
%         'LineStyle','-.','color','b','Marker','none','Linewidth',2);
% %     plot(linspace(min(VF_HotWT10uprimeshort(z,:)), max(VF_HotWT10uprimeshort(z,:)), 1000), pdf_u_fit, 'r-', 'LineWidth', 1.5); % Normal distribution
%     
%     legend(sprintf('Sub--window mean, SLPIV $\\mathrm{z}/\\delta$ = %.2f',VF_SLPIVASL.z(z,1)/VF_SLPIVASL.delta),...
%         sprintf('Global mean, SLPIV $\\mathrm{z}/\\delta$ = %.2f',VF_SLPIVASL.z(z,1)/VF_SLPIVASL.delta),...
%         'Interpreter','latex','FontSize',11,'Position',...
%         [0.236985524739487,0.810046174733969,0.556944479000946,0.096306066085292],...
%         'Numcolumns',1,'Orientation','vertical','color','none');
%     xlabel('$\mathrm{u}^{\prime}/u_{\tau}$','Interpreter','latex')
%     ylabel('p.d.f.','Interpreter','latex')
%     set(gca,'TickLabelInterpreter','latex','FontSize',14,'XGrid','on','YGrid','on')
%     xlim([-6.5 6.5])
%     ylim([0 0.3])
%     axis square
    
%     figure
%     plot(uprime_shorttime_ASL,pdf_uprime_shorttime_ASL,...
%         'LineStyle','none','color','r','Marker','^');
%     hold on
%     plot(linspace(min(VF_SLPIVASLuprimeshort(z,:,:),[],'all'), max(VF_SLPIVASLuprimeshort(z,:),[],'all'), 1000),...
%         pdf_u_fit, 'r-', 'LineWidth', 1.5); % Normal distribution
%     
%     figure
%     plot(wprime_shorttime_ASL,pdf_wprime_shorttime_ASL,...
%     'LineStyle','none','color','k','Marker','^');
%     hold on
%     plot(linspace(min(VF_SLPIVASLwprimeshort(z,:,:),[],'all'), max(VF_SLPIVASLwprimeshort(z,:,:),[],'all'), 1000),...
%         pdf_w_fit, 'k-', 'LineWidth', 1.5); % Normal distribution
    
end

figure 
plot(VF_HotWT7.z/VF_HotWT7.delta,std_uprime_shorttime_HotWT7,...
    'LineStyle','none','color','k','Marker','^')
hold on
plot(VF_HotWT10.z/VF_HotWT10.delta,std_uprime_shorttime_HotWT10,...
    'LineStyle','none','color','r','Marker','v')
plot(VF_SLPIVASL.z/VF_SLPIVASL.delta,std_uprime_shorttime_SLPIVASL,...
    'LineStyle','none','color','g','Marker','o')
legend('uprimeWT7','uprimeWT10','uprimeASL')
ylim([0 0.3])
figure 
plot(VF_HotWT7.z/VF_HotWT7.delta,std_wprime_shorttime_HotWT7,...
    'LineStyle','none','color','k','Marker','^')
hold on
plot(VF_HotWT10.z/VF_HotWT10.delta,std_wprime_shorttime_HotWT10,...
    'LineStyle','none','color','r','Marker','v')
plot(VF_SLPIVASL.z/VF_SLPIVASL.delta,std_wprime_shorttime_SLPIVASL,...
    'LineStyle','none','color','g','Marker','o')
legend('wprimeWT7','wprimeWT10','wprimeASL')
ylim([0 0.3])

%% Plot samples of velocity signal in the order of lambda_T



figure
plot(X_domain_WT7(1,:) / lambda_T_WT7(1,1),VF_HotWT7.u(1,:))
xlim([0 10])

figure
plot(T_domain_WT7(1,:),VF_HotWT7.u(1,:))

figure
plot((0:0.0001:1-0.0001).*mean(VF_HotWT7.u(1,:),2),VF_HotWT7.u(1,1:1:10000))
xlim([0.01 0.02])


figure
plot((0:1/120:1-1/120).*mean(VF_SLPIVASL.u(1,25,:),3),reshape(VF_SLPIVASL.u(1,25,1:1:120),[],1))

%% Interpolation to generated field

std_uprime_shorttime_GenWT7 = zeros(size(VF_GenWT7.z,1),1);
std_wprime_shorttime_GenWT7 = zeros(size(VF_GenWT7.z,1),1);
std_uprime_shorttime_GenWT10 = zeros(size(VF_GenWT10.z,1),1);
std_wprime_shorttime_GenWT10 = zeros(size(VF_GenWT10.z,1),1);
std_uprime_shorttime_GenASL = zeros(size(VF_GenASL.z,1),1);
std_wprime_shorttime_GenASL = zeros(size(VF_GenASL.z,1),1);

for z = 1:size(VF_GenWT7.z,1) 
    [r]=find(VF_GenWT7.z(z)<=VF_HotWT7.z,1,'first');
    std_uprime_shorttime_GenWT7(z) = std_uprime_shorttime_HotWT7(r);
    std_wprime_shorttime_GenWT7(z) = std_wprime_shorttime_HotWT7(r);
end

for z = 1:size(VF_GenWT10.z,1) 
    [r]=find(VF_GenWT10.z(z)<=VF_HotWT10.z,1,'first');
    std_uprime_shorttime_GenWT10(z) = std_uprime_shorttime_HotWT10(r);
    std_wprime_shorttime_GenWT10(z) = std_wprime_shorttime_HotWT10(r);
end

for z = 1:size(VF_GenASL.z,1) 
    [r]=find(VF_GenASL.z(z)<=VF_SLPIVASL.z,1,'first');
    std_uprime_shorttime_GenASL(z) = std_uprime_shorttime_SLPIVASL(r);
    std_wprime_shorttime_GenASL(z) = std_wprime_shorttime_SLPIVASL(r);
end
%% Saving
save('std_uprime_shorttime_Long_GenWT7.mat','std_uprime_shorttime_GenWT7')
save('std_uprime_shorttime_Long_GenWT10.mat','std_uprime_shorttime_GenWT10')
save('std_uprime_shorttime_Long_GenASL.mat','std_uprime_shorttime_GenASL')
save('std_wprime_shorttime_Long_GenWT7.mat','std_wprime_shorttime_GenWT7')
save('std_wprime_shorttime_Long_GenWT10.mat','std_wprime_shorttime_GenWT10')
save('std_wprime_shorttime_Long_GenASL.mat','std_wprime_shorttime_GenASL')