%% Data reading
Data3 = load('VF_HotWT7.mat');
VF_HotWT7 = Data3.VF_HotWT7;
clear('Data3')


Data6 = load('VF_HotWT10.mat');
VF_HotWT10 = Data6.VF_HotWT10;
clear('Data6')
%%
lambda_T_WT7 = zeros(size(VF_HotWT7.uprime, 1), 1);


for i = 1:size(VF_HotWT7.uprime, 1)
    
    current_row_epsilon = mean([VF_HotWT7.epsilon_str{i, :}{:}],2);
    
    lambda_T_WT7(i) = rms(VF_HotWT7.uprime(i,:),2)*sqrt(15*VF_HotWT7.nu/current_row_epsilon);
end


lambda_T_WT10 = zeros(size(VF_HotWT10.uprime, 1), 1);


for i = 1:size(VF_HotWT10.uprime, 1)
    
    current_row_epsilon = mean([VF_HotWT10.epsilon_str{i, :}{:}],2);
    
    lambda_T_WT10(i) = rms(VF_HotWT10.uprime(i,:),2)*sqrt(15*VF_HotWT10.nu/current_row_epsilon);
end

%% FFT

% Ns_HotWT7 = size(VF_HotWT7.u,2);
% lambda_T_time = lambda_T_WT7./mean(VF_HotWT7.u,2);
% n = ceil(Ns_HotWT7./(lambda_T_time.*VF_HotWT7.fs));
% 
% PSu_HotWT7 = fft(VF_HotWT7.u,[],2);
% PSu_short_HotWT7 = zeros(size(PSu_HotWT7));
% for z =1:size(n,1)
%     PSu_short_HotWT7(:,n(z,1):end) = PSu_HotWT7(:,n(z,1):end);
% end
% u_short_HotWT7=real(ifft(PSu_short_HotWT7,[],2));
% u_prime_short_HotWT7 = u_short_HotWT7;
% 
% 
% PSw_HotWT7 = fft(VF_HotWT7.w,[],2);
% PSw_short_HotWT7 = zeros(size(PSw_HotWT7));
% for z =1:size(n,1)
%     PSw_short_HotWT7(:,n(z,1):end) = PSw_HotWT7(:,n(z,1):end);
% end
% w_short_HotWT7=real(ifft(PSw_short_HotWT7,[],2));
% w_prime_short_HotWT7 = w_short_HotWT7 ;


Ns_HotWT10 = size(VF_HotWT10.u,2);
lambda_T_time = lambda_T_WT10./mean(VF_HotWT10.u,2);
n = ceil(Ns_HotWT10./(lambda_T_time.*VF_HotWT10.fs));

PSu_HotWT10 = fft(VF_HotWT10.u,[],2);
PSu_short_HotWT10 = zeros(size(PSu_HotWT10));
for z =1:size(n,1)
    PSu_short_HotWT10(:,n(z,1):end) = PSu_HotWT10(:,n(z,1):end);
end
u_short_HotWT10=real(ifft(PSu_short_HotWT10,[],2));
u_prime_short_HotWT10 = u_short_HotWT10;


PSw_HotWT10 = fft(VF_HotWT10.w,[],2);
PSw_short_HotWT10 = zeros(size(PSw_HotWT10));
for z =1:size(n,1)
    PSw_short_HotWT10(:,n(z,1):end) = PSw_HotWT10(:,n(z,1):end);
end
w_short_HotWT10=real(ifft(PSw_short_HotWT10,[],2));
w_prime_short_HotWT10 = w_short_HotWT10 ;

%% Analyzing

% mean_uprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
% std_uprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
% mean_wprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
% std_wprime_shorttime_HotWT7 = zeros(size(VF_HotWT7.z,1),1);
% 
% for z = 1: size(VF_HotWT7.z,1)
%     
%     
%     
%     [pdf_uprime_shorttime_WT7,uprime_shorttime_WT7] = ksdensity(u_prime_short_HotWT7(z,:),...
%     linspace(min(u_prime_short_HotWT7(z,:)), max(u_prime_short_HotWT7(z,:)), 50),'NumPoints',50);
% 
%     [pdf_wprime_shorttime_WT7,wprime_shorttime_WT7] = ksdensity(w_prime_short_HotWT7(z,:),...
%     linspace(min(w_prime_short_HotWT7(z,:)), max(w_prime_short_HotWT7(z,:)), 50),'NumPoints',50);
% 
% 
%     mean_uprime_shorttime_HotWT7(z,1) = mean(u_prime_short_HotWT7(z,:),2);
%     std_uprime_shorttime_HotWT7(z,1) = std(u_prime_short_HotWT7(z,:),0,2);
%     
%     mean_wprime_shorttime_HotWT7(z,1) = mean(w_prime_short_HotWT7(z,:),2);
%     std_wprime_shorttime_HotWT7(z,1) = std(w_prime_short_HotWT7(z,:),0,2);
%     
%     pdf_u_fit = normpdf(linspace(min(u_prime_short_HotWT7(z,:)), max(u_prime_short_HotWT7(z,:)), 1000), 0, std_uprime_shorttime_HotWT7(z,1));
%     pdf_w_fit = normpdf(linspace(min(w_prime_short_HotWT7(z,:)), max(w_prime_short_HotWT7(z,:)), 1000), 0, std_wprime_shorttime_HotWT7(z,1));
%     
%     figure
%     plot(uprime_shorttime_WT7,pdf_uprime_shorttime_WT7,...
%         'LineStyle','none','color','r','Marker','^');
%     hold on
%     plot(linspace(min(u_prime_short_HotWT7(z,:)), max(u_prime_short_HotWT7(z,:)), 1000), pdf_u_fit, 'r-', 'LineWidth', 1.5); % Normal distribution
%     
%     figure
%     plot(wprime_shorttime_WT7,pdf_wprime_shorttime_WT7,...
%     'LineStyle','none','color','k','Marker','^');
%     hold on
%     plot(linspace(min(w_prime_short_HotWT7(z,:)), max(w_prime_short_HotWT7(z,:)), 1000), pdf_w_fit, 'k-', 'LineWidth', 1.5); % Normal distribution
%     
% end



mean_uprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_uprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
mean_wprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);
std_wprime_shorttime_HotWT10 = zeros(size(VF_HotWT10.z,1),1);

for z = 1: size(VF_HotWT10.z,1)
    
    
    
    [pdf_uprime_shorttime_WT10,uprime_shorttime_WT10] = ksdensity(u_prime_short_HotWT10(z,:),...
    linspace(min(u_prime_short_HotWT10(z,:)), max(u_prime_short_HotWT10(z,:)), 50),'NumPoints',50);

    [pdf_wprime_shorttime_WT10,wprime_shorttime_WT10] = ksdensity(w_prime_short_HotWT10(z,:),...
    linspace(min(w_prime_short_HotWT10(z,:)), max(w_prime_short_HotWT10(z,:)), 50),'NumPoints',50);


    mean_uprime_shorttime_HotWT10(z,1) = mean(u_prime_short_HotWT10(z,:),2);
    std_uprime_shorttime_HotWT10(z,1) = std(u_prime_short_HotWT10(z,:),0,2);
    
    mean_wprime_shorttime_HotWT10(z,1) = mean(w_prime_short_HotWT10(z,:),2);
    std_wprime_shorttime_HotWT10(z,1) = std(w_prime_short_HotWT10(z,:),0,2);
    
    pdf_u_fit = normpdf(linspace(min(u_prime_short_HotWT10(z,:)), max(u_prime_short_HotWT10(z,:)), 1000), 0, std_uprime_shorttime_HotWT10(z,1));
    pdf_w_fit = normpdf(linspace(min(w_prime_short_HotWT10(z,:)), max(w_prime_short_HotWT10(z,:)), 1000), 0, std_wprime_shorttime_HotWT10(z,1));
    
    figure
    plot(uprime_shorttime_WT10,pdf_uprime_shorttime_WT10,...
        'LineStyle','none','color','r','Marker','^');
    hold on
    plot(linspace(min(u_prime_short_HotWT10(z,:)), max(u_prime_short_HotWT10(z,:)), 1000), pdf_u_fit, 'r-', 'LineWidth', 1.5); % Normal distribution
    
    figure
    plot(wprime_shorttime_WT10,pdf_wprime_shorttime_WT10,...
    'LineStyle','none','color','k','Marker','^');
    hold on
    plot(linspace(min(w_prime_short_HotWT10(z,:)), max(w_prime_short_HotWT10(z,:)), 1000), pdf_w_fit, 'k-', 'LineWidth', 1.5); % Normal distribution
    
end