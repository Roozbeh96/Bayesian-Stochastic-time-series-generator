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

mean_uprime_O_u_tau_chopped_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
std_uprime_O_u_tau_chopped_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
mean_wprime_O_u_tau_chopped_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
std_wprime_O_u_tau_chopped_HotWT7 = zeros(size(VF_HotWT7.z,1),1);

% Since the mean of each chopped signal is zero, there is no difference
% whether compute the std of each segment then take an average, or take std
% of the whole chopped signal. Here, I took std of the whole chopped
% signal.

for z = 1: size(VF_HotWT7.z,1)
    
    
    [pdf_uprime_chopped_WT7,uprime_chopped_WT7] = ksdensity(VF_HotWT7uprime_chopped{z}/VF_HotWT7.u_tau,...
    linspace(-6, 6, 100));

    [pdf_uprime_globaltime_WT7,uprime_globaltime_WT7] = ksdensity((VF_HotWT7.u(z,:)-mean(VF_HotWT7.u(z,:),2))/VF_HotWT7.u_tau,...
    linspace(-6, 6, 100));



    mean_uprime_O_u_tau_chopped_HotWT7(z,1) = mean(VF_HotWT7uprime_chopped{z}/VF_HotWT7.u_tau,2);
    std_uprime_O_u_tau_chopped_HotWT7(z,1) = std(VF_HotWT7uprime_chopped{z}/VF_HotWT7.u_tau,0,2);
    
    mean_wprime_O_u_tau_chopped_HotWT7(z,1) = mean(VF_HotWT7wprime_chopped{z}/VF_HotWT7.u_tau,2);
    std_wprime_O_u_tau_chopped_HotWT7(z,1) = std(VF_HotWT7wprime_chopped{z}/VF_HotWT7.u_tau,0,2);
    
    % Fig2(a) in the repo

    % figure
    % set(gcf,'Position',[622,508,806,394])
    % axes('Position',[0.08560794044665,0.134517766497462,0.40818858560794,0.83502538071066])
    % plot(uprime_chopped_WT7,pdf_uprime_chopped_WT7,...
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


mean_uprime_O_u_tau_chopped_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_uprime_O_u_tau_chopped_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
mean_wprime_O_u_tau_chopped_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_wprime_O_u_tau_chopped_HotWT10 = zeros(size(VF_HotWT10.z,1),1);

for z = 1: size(VF_HotWT10.z,1)
    
    
    
    [pdf_uprime_chopped_WT10,uprime_chopped_WT10] = ksdensity(VF_HotWT10uprime_chopped{z}/VF_HotWT10.u_tau,...
    linspace(-6, 6, 100));

    [pdf_uprime_globaltime_WT10,uprime_globaltime_WT10] = ksdensity((VF_HotWT10.u(z,:)-mean(VF_HotWT10.u(z,:),2))/VF_HotWT10.u_tau,...
    linspace(-6, 6, 100));


    mean_uprime_O_u_tau_chopped_HotWT10(z,1) = mean(VF_HotWT10uprime_chopped{z}/VF_HotWT10.u_tau,2);
    std_uprime_O_u_tau_chopped_HotWT10(z,1) = std(VF_HotWT10uprime_chopped{z}/VF_HotWT10.u_tau,0,2);
    
    mean_wprime_O_u_tau_chopped_HotWT10(z,1) = mean(VF_HotWT10wprime_chopped{z}/VF_HotWT10.u_tau,2);
    std_wprime_O_u_tau_chopped_HotWT10(z,1) = std(VF_HotWT10wprime_chopped{z}/VF_HotWT10.u_tau,0,2);
    

    %Fig2(b) in the repo

    % figure
    % set(gcf,'Position',[622,508,806,394])
    % axes('Position',[0.08560794044665,0.134517766497462,0.40818858560794,0.83502538071066])
    % plot(uprime_chopped_WT10,pdf_uprime_chopped_WT10,...
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
%% Hypothesis Testing. 

%Since the statistics of HotWT7 and HotWT10 are similar to each other, we just use
% the statistics of HotWT7 for signal generation. But before that, we need
% to justify our claim. To do so, I use hypothesis testing.
%H_0: std of chopped velocity signal at z = 5[cm](ind = 9) for both HotWT7
%and HotWT10 have difference MORE than 0.01(margin)-> HotWT7 and HotWT10
%are NOT equivalent->(NOT FAVORABLE)->|\sigma_{WT7}-\sigma_{WT10}| > margin.
%H_a: std of chopped velocity signal at z = 5[cm](ind = 9) for both HotWT7
%and HotWT10 have difference LESS than 0.01(margin)-> HotWT7 and HotWT10
%are equivalent->(FAVORABLE)->|\sigma_{WT7}-\sigma_{WT10}| <= margin.
%Always put FAVORABLE option in H_a.
% Since we do not know about the distribution (Gaussian or other ones), I
% used non-parametric method for hypothesis testing. Use bootstraping
% technique for sampling. 

margin = 0.025;
x = VF_HotWT7uprime_chopped{9}/VF_HotWT7.u_tau;
y = VF_HotWT10uprime_chopped{9}/VF_HotWT10.u_tau;
n = numel(x); m = numel(y);


% Now, let us bootstraping to from x, and y to generate samples of xb and 
% yb and compute the difference between the stds (d_boot). Having d_boot,
% we can have the distribution of the standard deviation difference. Using
% this distribution, we can find the confidence interval as well. I 
% repeated this scenario 20,000 times. Distribution of d_boot is our 
% criterion because it is calculated based on samples which is organic.


B = 20000;
d_boot = zeros(B,1);
for b = 1:B
    xb = x(randi(n, floor(n/10), 1));
    yb = y(randi(m, floor(m/10), 1));
    d_boot(b) = std(xb,0,2) - std(yb,0,2);
end

alpha = 0.05; % significance level

% Percentile CI for d_boot at level (1-2*alpha). This means that, based on
% the distribution, with the probability of (1-2*alpha)*100 percent, the
% difference between standard deviation is between [lo hi](BASED ON 
% DISTRIBUTION). As the alpha increases, the confidence interval reduces
% and the probability of rejecting H0 increases. if CI be in the range of
% (-margin,margin)->Favorable case-> we can reject the null hypothesis. 
% If 'reject_H0 = 1' means we can reject the null-hypothesis.

lo = prctile(d_boot, 100*alpha);
hi = prctile(d_boot, 100*(1-alpha));
fprintf('reject_H0: %d\n',lo>=-margin & hi<=margin)

% Other approach than we can test our hypothesis, is using p-value.
% Calculate the p-value based on the assumption of NULL HYPOTHESIS(H_0). If
% p-val be less than alpha, which shows the probability of occurrence of
% null hypothesis based on the distribution, which is always our criteria,
% is low-> So we can reject the H_0.

pval_lo = mean(d_boot < -margin);
pval_hi = mean(d_boot > margin);

fprintf('reject_H0: %d\n',pval_lo < alpha & pval_hi< alpha)

%% Interpolation to generated field


z_min = 50*VF_HotWT7.nu / VF_HotWT7.u_tau;
z_max = 0.25*VF_HotWT7.delta;
lambda_T = 0.01;

VF_GenWT7_z = z_min:0.4*lambda_T/10:z_max + 0.4*lambda_T/10;


std_uprime_O_u_tau_chopped = zeros(size(VF_GenWT7_z,2),1);
std_wprime_O_u_tau_chopped = zeros(size(VF_GenWT7_z,2),1);


for z = 1:size(VF_GenWT7_z,2) 
    [r]=find(VF_GenWT7_z(z)<=VF_HotWT7.z,1,'first');
    std_uprime_O_u_tau_chopped(z) = std_uprime_O_u_tau_chopped_HotWT7(r);
    std_wprime_O_u_tau_chopped(z) = std_wprime_O_u_tau_chopped_HotWT7(r);
end


%% Saving

save('std_uprime_O_u_tau_chopped.mat','std_uprime_O_u_tau_chopped')
save('std_wprime_O_u_tau_chopped.mat','std_wprime_O_u_tau_chopped')

