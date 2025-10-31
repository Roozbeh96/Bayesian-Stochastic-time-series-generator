%#ok<*NOPTS>

%{
    In this code, the experimental velocity signal is chopped into desired
    length-scale.
%}
%% Read data.

Data = load('VF_HotWT7.mat');
VF_HotWT7 = Data.VF_HotWT7;
clear('Data')


Data = load('VF_HotWT10.mat');
VF_HotWT10 = Data.VF_HotWT10;
clear('Data')


%% Chopping the signal into L

T_domain_WT7 = (0:1:size(VF_HotWT7.u,2)-1)*1/VF_HotWT7.fs;
X_domain_WT7 = mean(VF_HotWT7.u,2)*T_domain_WT7;

T_domain_WT10 = (0:1:size(VF_HotWT10.u,2)-1)*1/VF_HotWT10.fs;
X_domain_WT10 = mean(VF_HotWT10.u,2)*T_domain_WT10;

VF_HotWT7uprime_chopped = cell(size(VF_HotWT7.z,1),1);
VF_HotWT7wprime_chopped = cell(size(VF_HotWT7.z,1),1);
VF_HotWT10uprime_chopped = cell(size(VF_HotWT10.z,1),1);
VF_HotWT10wprime_chopped = cell(size(VF_HotWT10.z,1),1);


LWT7 = 0.7;

for z = 1: size(VF_HotWT7.z,1)
    

    ratio_WT7= X_domain_WT7(z,:) / LWT7;
    unitindex_WT7 = find(ratio_WT7 >= 1, 1, 'first');
    unitindex_WT7 = unitindex_WT7 - 1;
    
    for s = 1 : floor(size(VF_HotWT7.u(z,:),2)/unitindex_WT7)

        VF_HotWT7uprime_chopped{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2);

        VF_HotWT7wprime_chopped{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2);
        
        %{
            Run the code below if you want Fig 1 in the repo.
        %}

        % plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/LWT7,...
        %     VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
        % hold on
        % plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/LWT7,...
        %     VF_HotWT7uprime_chopped{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
        % 
        % plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/LWT7,...
        %     VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-mean(VF_HotWT7.u(z,:),2),'linewidth',2)
        % legend(sprintf('u[m/s],H-W(m1) z$/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
        %     sprintf('$\\mathrm{u}^{\\prime}_{\\mathrm{Local}}\\mathrm{[m/s]}$,H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
        %     sprintf('$\\mathrm{u}^{\\prime}_{\\mathrm{Global}}\\mathrm{[m/s]}$,H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
        %     'Interpreter','latex','FontSize',10,'Position',...
        %     [0.69418761794644,0.810046174733969,0.271111721158469,0.096306066085292],...
        %     'Numcolumns',1,'Orientation','vertical','color','none');
        % xlabel('x/L','Interpreter','latex')
        % ylabel('$\mathrm{u},\mathrm{u}^{\prime}$[m/s]','Interpreter','latex')
        % set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on','YGrid','on')

    end

    
end


for z=1:size(VF_HotWT7.z,1)
    fprintf('At z/\\delta = %.2f\n',VF_HotWT7.z(z)/VF_HotWT7.delta)
    fprintf('Corrcoef(u,w)= %.2f\n',corrcoef(VF_HotWT7uprime_chopped{z},VF_HotWT7wprime_chopped{z}))
    
end


LWT10 = 0.7;
for z = 1: size(VF_HotWT10.z,1)
    

    ratio_WT10 = X_domain_WT10(z,:) / LWT10;
    unitindex_WT10 = find(ratio_WT10 >= 1, 1, 'first');
    unitindex_WT10 = unitindex_WT10 - 1;
    
    for s = 1 : floor(size(VF_HotWT10.u(z,:),2)/unitindex_WT10)
        VF_HotWT10uprime_chopped{z}((s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2);
        
        VF_HotWT10wprime_chopped{z}((s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2);        
    end
    
end


for z=1:size(VF_HotWT10.z,1)
    
    fprintf('At z/\\delta = %.2f\n',VF_HotWT10.z(z)/VF_HotWT10.delta)
    fprintf('Corrcoef(u,w)= %.2f\n',corrcoef(VF_HotWT10uprime_chopped{z},VF_HotWT10wprime_chopped{z}))
    
end


%% Analyzing the distriburtion

mean_uprime_O_u_tau_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
std_uprime_O_u_tau_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
mean_wprime_O_u_tau_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
std_wprime__O_u_tau_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);


for z = 1: size(VF_HotWT7.z,1)
    
    
    [pdf_uprime_shorttime_WT7,uprime_shorttime_WT7] = ksdensity(VF_HotWT7uprime_chopped{z}/VF_HotWT7.u_tau,...
    linspace(-6, 6, 100));

    [pdf_uprime_globaltime_WT7,uprime_globaltime_WT7] = ksdensity((VF_HotWT7.u(z,:)-mean(VF_HotWT7.u(z,:),2))/VF_HotWT7.u_tau,...
    linspace(-6, 6, 100));



    mean_uprime_O_u_tau_shorttime_HotWT7(z,1) = mean(VF_HotWT7uprime_chopped{z}/VF_HotWT7.u_tau,2);
    std_uprime_O_u_tau_shorttime_HotWT7(z,1) = std(VF_HotWT7uprime_chopped{z}/VF_HotWT7.u_tau,0,2);
    
    mean_wprime_O_u_tau_shorttime_HotWT7(z,1) = mean(VF_HotWT7wprime_chopped{z}/VF_HotWT7.u_tau,2);
    std_wprime__O_u_tau_shorttime_HotWT7(z,1) = std(VF_HotWT7wprime_chopped{z}/VF_HotWT7.u_tau,0,2);
    
    % Fig2(a) in the repo

    % figure
    % set(gcf,'Position',[622,508,806,394])
    % axes('Position',[0.08560794044665,0.134517766497462,0.40818858560794,0.83502538071066])
    % plot(uprime_shorttime_WT7,pdf_uprime_shorttime_WT7,...
    %     'LineStyle','-','color','r','Marker','none','Linewidth',2);
    % hold on
    % plot(uprime_globaltime_WT7,pdf_uprime_globaltime_WT7,...
    %     'LineStyle','-.','color','b','Marker','none','Linewidth',2);
    % 
    % legend(sprintf('Sub--window mean, H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
    %     sprintf('Global mean, H-W(m1) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT7.z(z,1)/VF_HotWT7.delta),...
    %     'Interpreter','latex','FontSize',11,'Position',...
    %     [0.094243393010764,0.853855920868091,0.389991769448501,0.100761418173156],...
    %     'Numcolumns',1,'Orientation','vertical','color',[1,1,1]);
    % xlabel('$\mathrm{u}^{\prime}/u_{\tau}$','Interpreter','latex')
    % ylabel('p.d.f.','Interpreter','latex')
    % annotation(gcf,'textbox',...
    % [-0.00893788819875776 0.926395939086293 0.0412360248447205 0.088832487309642],...
    % 'String','(a)',...
    % 'Interpreter','latex',...
    % 'FontSize',18,...
    % 'FitBoxToText','off',...
    % 'EdgeColor','none');
    % set(gca,'TickLabelInterpreter','latex','FontSize',14,'XGrid','on','YGrid','on')
    % xlim([-5.5 5.5])
    % ylim([0 0.3])
    % axis square
    

end


mean_uprime_O_u_tau_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_uprime_O_u_tau_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
mean_wprime_O_u_tau_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_wprime_O_u_tau_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);

for z = 1: size(VF_HotWT10.z,1)
    
    
    
    [pdf_uprime_shorttime_WT10,uprime_shorttime_WT10] = ksdensity(VF_HotWT10uprime_chopped{z}/VF_HotWT10.u_tau,...
    linspace(-6, 6, 100));

    [pdf_uprime_globaltime_WT10,uprime_globaltime_WT10] = ksdensity((VF_HotWT10.u(z,:)-mean(VF_HotWT10.u(z,:),2))/VF_HotWT10.u_tau,...
    linspace(-6, 6, 100));


    mean_uprime_O_u_tau_shorttime_HotWT10(z,1) = mean(VF_HotWT10uprime_chopped{z}/VF_HotWT10.u_tau,2);
    std_uprime_O_u_tau_shorttime_HotWT10(z,1) = std(VF_HotWT10uprime_chopped{z}/VF_HotWT10.u_tau,0,2);
    
    mean_wprime_O_u_tau_shorttime_HotWT10(z,1) = mean(VF_HotWT10wprime_chopped{z}/VF_HotWT10.u_tau,2);
    std_wprime_O_u_tau_shorttime_HotWT10(z,1) = std(VF_HotWT10wprime_chopped{z}/VF_HotWT10.u_tau,0,2);
    

    %Fig2(b) in the repo
    % figure
    % set(gcf,'Position',[622,508,806,394])
    % axes('Position',[0.08560794044665,0.134517766497462,0.40818858560794,0.83502538071066])
    % plot(uprime_shorttime_WT10,pdf_uprime_shorttime_WT10,...
    %     'LineStyle','-','color','r','Marker','none','Linewidth',2);
    % hold on
    % plot(uprime_globaltime_WT10,pdf_uprime_globaltime_WT10,...
    %     'LineStyle','-.','color','b','Marker','none','Linewidth',2);
    % legend(sprintf('Sub--window mean, H-W(m2) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT10.z(z,1)/VF_HotWT10.delta),...
    %     sprintf('Global mean, H-W(m2) $\\mathrm{z}/\\delta$ = %.2f',VF_HotWT10.z(z,1)/VF_HotWT10.delta),...
    %     'Interpreter','latex','FontSize',11,'Position',...
    %     [0.094243393010764,0.853855920868091,0.389991769448501,0.100761418173156],...
    %     'Numcolumns',1,'Orientation','vertical','color',[1,1,1]);
    % xlabel('$\mathrm{u}^{\prime}/u_{\tau}$','Interpreter','latex')
    % ylabel('p.d.f.','Interpreter','latex')
    % annotation(gcf,'textbox',...
    % [0.486714285714286 0.918781725888321 0.0412360248447204 0.088832487309642],...
    % 'String','(b)',...
    % 'Interpreter','latex',...
    % 'FontSize',18,...
    % 'FitBoxToText','off',...
    % 'EdgeColor','none');
    % set(gca,'TickLabelInterpreter','latex','FontSize',14,'XGrid','on','YGrid','on')
    % xlim([-5.5 5.5])
    % ylim([0 0.3])
    % axis square


end



figure 
plot(VF_HotWT7.z/VF_HotWT7.delta,std_uprime_O_u_tau_shorttime_HotWT7,...
    'LineStyle','none','color','k','Marker','^')
hold on
plot(VF_HotWT10.z/VF_HotWT10.delta,std_uprime_O_u_tau_shorttime_HotWT10,...
    'LineStyle','none','color','r','Marker','v')
legend('uprimeWT7','uprimeWT10')
ylim([0 0.3])
figure 
plot(VF_HotWT7.z/VF_HotWT7.delta,std_wprime__O_u_tau_shorttime_HotWT7,...
    'LineStyle','none','color','k','Marker','^')
hold on
plot(VF_HotWT10.z/VF_HotWT10.delta,std_wprime_O_u_tau_shorttime_HotWT10,...
    'LineStyle','none','color','r','Marker','v')
legend('wprimeWT7','wprimeWT10')
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
    std_uprime_shorttime_GenWT7(z) = std_uprime_O_u_tau_shorttime_HotWT7(r);
    std_wprime_shorttime_GenWT7(z) = std_wprime__O_u_tau_shorttime_HotWT7(r);
end

for z = 1:size(VF_GenWT10.z,1) 
    [r]=find(VF_GenWT10.z(z)<=VF_HotWT10.z,1,'first');
    std_uprime_shorttime_GenWT10(z) = std_uprime_O_u_tau_shorttime_HotWT10(r);
    std_wprime_shorttime_GenWT10(z) = std_wprime_O_u_tau_shorttime_HotWT10(r);
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