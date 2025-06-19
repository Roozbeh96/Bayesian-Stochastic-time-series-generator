Data3 = load('VF_HotWT7.mat');
VF_HotWT7 = Data3.VF_HotWT7;
clear('Data3')


Data6 = load('VF_HotWT10.mat');
VF_HotWT10 = Data6.VF_HotWT10;
clear('Data6')


Data8 = load('VF_SLPIVASL.mat');
VF_SLPIVASL = Data8.VF_SLPIVASL;
clear('Data8')
%%

T_domain_WT7 = (0:1:size(VF_HotWT7.u,2)-1)*1/VF_HotWT7.fs;
X_domain_WT7 = mean(VF_HotWT7.u,2)*T_domain_WT7;
res_WT7 = X_domain_WT7(:,2)-X_domain_WT7(:,1);

T_domain_WT10 = (0:1:size(VF_HotWT10.u,2)-1)*1/VF_HotWT10.fs;
X_domain_WT10 = mean(VF_HotWT10.u,2)*T_domain_WT10;
res_WT10 = X_domain_WT10(:,2)-X_domain_WT10(:,1);

T_domain_ASL = (0:1:size(VF_SLPIVASL.u,3)-1)*1/VF_SLPIVASL.fs;
T_domain_ASL = reshape(T_domain_ASL,1,1,[]);
X_domain_ASL = mean(VF_SLPIVASL.u,3).*T_domain_ASL;
res_ASL = X_domain_ASL(:,:,2)-X_domain_ASL(:,:,1);

VF_HotWT7uprimelong = VF_HotWT7.u - mean(VF_HotWT7.u,2);
VF_HotWT7wprimelong = VF_HotWT7.w;
VF_HotWT10uprimelong = VF_HotWT10.u - mean(VF_HotWT10.u,2);
VF_HotWT10wprimelong = VF_HotWT10.w;
VF_SLPIVASLuprimelong = VF_SLPIVASL.uprime;
VF_SLPIVASLwprimelong = VF_SLPIVASL.w;


%% Choping the signal into lambda_T

T_domain_WT7 = (0:1:size(VF_HotWT7.u,2)-1)*1/VF_HotWT7.fs;
X_domain_WT7 = mean(VF_HotWT7.u,2)*T_domain_WT7;
res_WT7 = X_domain_WT7(:,2)-X_domain_WT7(:,1);

T_domain_WT10 = (0:1:size(VF_HotWT10.u,2)-1)*1/VF_HotWT10.fs;
X_domain_WT10 = mean(VF_HotWT10.u,2)*T_domain_WT10;
res_WT10 = X_domain_WT10(:,2)-X_domain_WT10(:,1);

T_domain_ASL = (0:1:size(VF_SLPIVASL.u,3)-1)*1/VF_SLPIVASL.fs;
T_domain_ASL = reshape(T_domain_ASL,1,1,[]);
X_domain_ASL = mean(VF_SLPIVASL.u,3).*T_domain_ASL;
res_ASL = X_domain_ASL(:,:,2)-X_domain_ASL(:,:,1);

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



for z = 1: size(VF_HotWT7.z,1)
    
    %     ratio_WT7= X_domain_WT7(z,:) / lambda_T_WT7(z,1);
    ratio_WT7= X_domain_WT7(z,:) / 0.3;
    unitindex_WT7 = find(ratio_WT7 >= 1, 1, 'first');
    unitindex_WT7 = unitindex_WT7 - 1;
    
    for s = 1 : floor(size(VF_HotWT7.u(z,:),2)/unitindex_WT7)
        VF_HotWT7uprimeshort(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2);
        
        VF_HotWT7wprimeshort(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2);
        
        
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

for z = 1: size(VF_HotWT10.z,1)
    
%     ratio_WT10 = X_domain_WT10(z,:) / lambda_T_WT10(z,1);
    ratio_WT10 = X_domain_WT10(z,:) / 0.3;
    unitindex_WT10 = find(ratio_WT10 >= 1, 1, 'first');
    unitindex_WT10 = unitindex_WT10 - 1;
    
    for s = 1 : floor(size(VF_HotWT10.u(z,:),2)/unitindex_WT10)
        VF_HotWT10uprimeshort(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2);
        
        VF_HotWT10wprimeshort(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2);        
    end
    
end

cols_with_zeros_uprime_WT10 = any(VF_HotWT10uprimeshort == 0, 1);
VF_HotWT10uprimeshort = VF_HotWT10uprimeshort(:, ~cols_with_zeros_uprime_WT10);

cols_with_zeros_wprime_WT10 = any(VF_HotWT10wprimeshort == 0, 1);
VF_HotWT10wprimeshort = VF_HotWT10wprimeshort(:, ~cols_with_zeros_wprime_WT10);

%% min-max distribution

value_pairs_cell_HotWT7 = cell(size(VF_HotWT7.z,1), 1);

for z = 1:size(VF_HotWT7.z,1)
    % Find local maxima and minima for the current elevation
    [peaks, peak_locs] = findpeaks(VF_HotWT7uprimelong(z,:)); % Local maxima from long should be changed to short
    [inverted_peaks, min_locs] = findpeaks(-VF_HotWT7uprimelong(z,:)); % Local minima
    valleys = -inverted_peaks;

    % Initialize an empty matrix for storing value pairs at this elevation
    value_pairs = [];

    % Check if there are any maxima or minima found
    if isempty(peak_locs) || isempty(min_locs)
        disp(['No local maxima or minima found at elevation z = ', num2str(z)]);
    else
        % Start from the first max or min
        i = 1; % index for max_locs
        j = 1; % index for min_locs
        
        % Determine the correct sequence
        if peak_locs(1) < min_locs(1)
            % Case: First local maximum comes before first local minimum
            while i <= length(peak_locs) && j <= length(min_locs)
                value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of min]
                i = i + 1;
                if i <= length(peak_locs) % Only proceed if i is valid
                    value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of next max, value of min]
                    j = j + 1;
                end
            end
        else
            % Case: First local minimum comes before first local maximum
            while i <= length(peak_locs) && j <= length(min_locs)
                value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of min]
                j = j + 1;
                if j <= length(min_locs) % Only proceed if j is valid
                    value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of next min]
                    i = i + 1;
                end
            end
        end
    end

    % Store the value pairs in the cell array for this elevation
    value_pairs_cell_HotWT7{z} = value_pairs;


end



value_pairs_cell_HotWT10 = cell(size(VF_HotWT10.z,1), 1);

for z = 1:size(VF_HotWT10.z,1)
    % Find local maxima and minima for the current elevation
    [peaks, peak_locs] = findpeaks(VF_HotWT10uprimelong(z,:)); % Local maxima
    [inverted_peaks, min_locs] = findpeaks(-VF_HotWT10uprimelong(z,:)); % Local minima
    valleys = -inverted_peaks;

    % Initialize an empty matrix for storing value pairs at this elevation
    value_pairs = [];

    % Check if there are any maxima or minima found
    if isempty(peak_locs) || isempty(min_locs)
        disp(['No local maxima or minima found at elevation z = ', num2str(z)]);
    else
        % Start from the first max or min
        i = 1; % index for max_locs
        j = 1; % index for min_locs
        
        % Determine the correct sequence
        if peak_locs(1) < min_locs(1)
            % Case: First local maximum comes before first local minimum
            while i <= length(peak_locs) && j <= length(min_locs)
                value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of min]
                i = i + 1;
                if i <= length(peak_locs) % Only proceed if i is valid
                    value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of next max, value of min]
                    j = j + 1;
                end
            end
        else
            % Case: First local minimum comes before first local maximum
            while i <= length(peak_locs) && j <= length(min_locs)
                value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of min]
                j = j + 1;
                if j <= length(min_locs) % Only proceed if j is valid
                    value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of next min]
                    i = i + 1;
                end
            end
        end
    end

    % Store the value pairs in the cell array for this elevation
    value_pairs_cell_HotWT10{z} = value_pairs;


end



value_pairs_cell_SLPIVASL = cell(size(VF_SLPIVASL.z,1), 1);

for z = 1:size(VF_SLPIVASL.z,1)
    % Find local maxima and minima for the current elevation
    [peaks, peak_locs] = findpeaks(reshape(VF_SLPIVASLuprimelong(z,25,:),1,[])); % Local maxima
    [inverted_peaks, min_locs] = findpeaks(reshape(-VF_SLPIVASLuprimelong(z,25,:),1,[])); % Local minima
    valleys = -inverted_peaks;

    % Initialize an empty matrix for storing value pairs at this elevation
    value_pairs = [];

    % Check if there are any maxima or minima found
    if isempty(peak_locs) || isempty(min_locs)
        disp(['No local maxima or minima found at elevation z = ', num2str(z)]);
    else
        % Start from the first max or min
        i = 1; % index for max_locs
        j = 1; % index for min_locs
        
        % Determine the correct sequence
        if peak_locs(1) < min_locs(1)
            % Case: First local maximum comes before first local minimum
            while i <= length(peak_locs) && j <= length(min_locs)
                value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of min]
                i = i + 1;
                if i <= length(peak_locs) % Only proceed if i is valid
                    value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of next max, value of min]
                    j = j + 1;
                end
            end
        else
            % Case: First local minimum comes before first local maximum
            while i <= length(peak_locs) && j <= length(min_locs)
                value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of min]
                j = j + 1;
                if j <= length(min_locs) % Only proceed if j is valid
                    value_pairs = [value_pairs; peaks(i), valleys(j)]; % [value of max, value of next min]
                    i = i + 1;
                end
            end
        end
    end

    % Store the value pairs in the cell array for this elevation
    value_pairs_cell_SLPIVASL{z} = value_pairs;


end

%% Saving & Loading
save('value_pairs_cell_HotWT7.mat','value_pairs_cell_HotWT7')
save('value_pairs_cell_HotWT10.mat','value_pairs_cell_HotWT10')
save('value_pairs_cell_SLPIVASL.mat','value_pairs_cell_SLPIVASL')

load('value_pairs_cell_HotWT7.mat')
load('value_pairs_cell_HotWT10.mat')
load('value_pairs_cell_SLPIVASL.mat')
%% ploting

nbins = 1000;
data_HotWT7_1 = [value_pairs_cell_HotWT7{1}(:,2)/VF_HotWT7.u_tau, value_pairs_cell_HotWT7{1}(:,1)/VF_HotWT7.u_tau]; % Example correlated non-Gaussian data

% Compute 2D histogram (counts per bin)
[counts, X_edges, Y_edges] = histcounts2(data_HotWT7_1(:,1), data_HotWT7_1(:,2), nbins);

% Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = (Y_edges(1:end-1) + Y_edges(2:end)) / 2;
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (assumes uniform bins)
dx = mean(diff(X_edges)); % Bin width in X
dy = mean(diff(Y_edges)); % Bin width in Y
bin_area = dx * dy; % Area of each bin

% Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT7_1 = counts / (sum(counts(:)) * bin_area);

% Compute the total probability mass (should be ≈ 1)
total_probability = sum(pdf_hist_HotWT7_1(:)) * bin_area;
disp(['Total Probability Mass (should be ≈ 1): ', num2str(total_probability)]);

figure;

set(gcf,'Position',[800,317,797,582])
axes('Position',[0.066478460895022,0.506680943904672,0.278565453785027,0.41875])
scatter(data_HotWT7_1(:,1), data_HotWT7_1(:,2), 5, [0.0,1,0.0],'filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT7_1', 100, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 0.5])
hcb2 = colorbar;
hcb2.Location = "north";
hcb2.Position = [0.095357590966123,0.924460431654677,0.885821831869511,0.025903057081893];
title(hcb2, 'p.d.f.', 'Interpreter', 'Latex', 'FontSize', 20, 'position',...
    [-29.625000000000938,-4.775,0])
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation = 'out';
axis equal
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
% xlabel('Local min u$^{\prime}$','Interpreter','latex','FontSize',15)
ylabel('Local max u$^{\prime}/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
text(0.018018018018018,0.932432432432432,0,'(a)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
xlim([-6 6])
ylim([-6 6])

data_HotWT7_9 = [value_pairs_cell_HotWT7{9}(:,2)/VF_HotWT7.u_tau, value_pairs_cell_HotWT7{9}(:,1)/VF_HotWT7.u_tau]; % Example correlated non-Gaussian data
[counts, X_edges, Y_edges] = histcounts2(data_HotWT7_9(:,1), data_HotWT7_9(:,2), nbins);

% Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = (Y_edges(1:end-1) + Y_edges(2:end)) / 2;
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (assumes uniform bins)
dx = mean(diff(X_edges)); % Bin width in X
dy = mean(diff(Y_edges)); % Bin width in Y
bin_area = dx * dy; % Area of each bin

% Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT7_9 = counts / (sum(counts(:)) * bin_area);

% Compute the total probability mass (should be ≈ 1)
total_probability = sum(pdf_hist_HotWT7_9(:)) * bin_area;
disp(['Total Probability Mass (should be ≈ 1): ', num2str(total_probability)]);


axes('Position',[0.066478460895023,0.072022485104706,0.278565453785027,0.418750000000002])
scatter(data_HotWT7_9(:,1), data_HotWT7_9(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT7_9', 100, 'LineWidth', 3); % Contour of PDF
% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 0.5])
axis equal
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
text(0.018018018018018,0.932432432432432,0,'(b)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
xlabel('Local min u$^{\prime}/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
ylabel('Local max u$^{\prime}/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
xlim([-6 6])
ylim([-6 6])


data_HotWT10_2 = [value_pairs_cell_HotWT10{2}(:,2)/VF_HotWT10.u_tau, value_pairs_cell_HotWT10{2}(:,1)/VF_HotWT10.u_tau]; % Example correlated non-Gaussian data
[counts, X_edges, Y_edges] = histcounts2(data_HotWT10_2(:,1), data_HotWT10_2(:,2), nbins);

% Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = (Y_edges(1:end-1) + Y_edges(2:end)) / 2;
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (assumes uniform bins)
dx = mean(diff(X_edges)); % Bin width in X
dy = mean(diff(Y_edges)); % Bin width in Y
bin_area = dx * dy; % Area of each bin

% Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT10_2 = counts / (sum(counts(:)) * bin_area);

% Compute the total probability mass (should be ≈ 1)
total_probability = sum(pdf_hist_HotWT10_2(:)) * bin_area;
disp(['Total Probability Mass (should be ≈ 1): ', num2str(total_probability)]);


axes('Position',[0.38266415725637,0.506680943904679,0.278565453785019,0.418750000000001])
scatter(data_HotWT10_2(:,1), data_HotWT10_2(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT10_2', 100, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 0.5])
axis equal
text(0.018018018018018,0.932432432432432,0,'(c)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
xlim([-6 6])
ylim([-6 6])

data_HotWT10_9 = [value_pairs_cell_HotWT10{9}(:,2)/VF_HotWT10.u_tau, value_pairs_cell_HotWT10{9}(:,1)/VF_HotWT10.u_tau]; % Example correlated non-Gaussian data
[counts, X_edges, Y_edges] = histcounts2(data_HotWT10_9(:,1), data_HotWT10_9(:,2), nbins);

% Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = (Y_edges(1:end-1) + Y_edges(2:end)) / 2;
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (assumes uniform bins)
dx = mean(diff(X_edges)); % Bin width in X
dy = mean(diff(Y_edges)); % Bin width in Y
bin_area = dx * dy; % Area of each bin

% Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT10_9 = counts / (sum(counts(:)) * bin_area);

% Compute the total probability mass (should be ≈ 1)
total_probability = sum(pdf_hist_HotWT10_9(:)) * bin_area;
disp(['Total Probability Mass (should be ≈ 1): ', num2str(total_probability)]);


axes('Position',[0.38266415725637,0.071799982694269,0.278565453785019,0.418750000000002])
scatter(data_HotWT10_9(:,1), data_HotWT10_9(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT10_9', 100, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 0.5])
axis equal
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
text(0.018018018018018,0.932432432432432,0,'(d)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
xlabel('Local min u$^{\prime}/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
xlim([-6 6])
ylim([-6 6])

data_SLPIVASL_7 = [value_pairs_cell_SLPIVASL{7}(:,2)/VF_SLPIVASL.u_tau, value_pairs_cell_SLPIVASL{7}(:,1)/VF_SLPIVASL.u_tau]; % Example correlated non-Gaussian data
[counts, X_edges, Y_edges] = histcounts2(data_SLPIVASL_7(:,1), data_SLPIVASL_7(:,2), nbins);

% Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = (Y_edges(1:end-1) + Y_edges(2:end)) / 2;
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (assumes uniform bins)
dx = mean(diff(X_edges)); % Bin width in X
dy = mean(diff(Y_edges)); % Bin width in Y
bin_area = dx * dy; % Area of each bin

% Normalize histogram to create probability density function (PDF)
pdf_hist_SLPIVASL_7 = counts / (sum(counts(:)) * bin_area);

% Compute the total probability mass (should be ≈ 1)
total_probability = sum(pdf_hist_SLPIVASL_7(:)) * bin_area;
disp(['Total Probability Mass (should be ≈ 1): ', num2str(total_probability)]);


axes('Position',[0.703868674194876,0.506680943904679,0.278565453785017,0.418750000000001])
scatter(data_SLPIVASL_7(:,1), data_SLPIVASL_7(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_SLPIVASL_7', 100, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 0.5])
axis equal
text(0.018018018018018,0.932432432432432,0,'(e)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
xlim([-6 6])
ylim([-6 6])


data_SLPIVASL_19 = [value_pairs_cell_SLPIVASL{19}(:,2)/VF_SLPIVASL.u_tau, value_pairs_cell_SLPIVASL{19}(:,1)/VF_SLPIVASL.u_tau]; % Example correlated non-Gaussian data
[counts, X_edges, Y_edges] = histcounts2(data_SLPIVASL_19(:,1), data_SLPIVASL_19(:,2), nbins);

% Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = (Y_edges(1:end-1) + Y_edges(2:end)) / 2;
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (assumes uniform bins)
dx = mean(diff(X_edges)); % Bin width in X
dy = mean(diff(Y_edges)); % Bin width in Y
bin_area = dx * dy; % Area of each bin

% Normalize histogram to create probability density function (PDF)
pdf_hist_SLPIVASL_19 = counts / (sum(counts(:)) * bin_area);

% Compute the total probability mass (should be ≈ 1)
total_probability = sum(pdf_hist_SLPIVASL_19(:)) * bin_area;
disp(['Total Probability Mass (should be ≈ 1): ', num2str(total_probability)]);


axes('Position',[0.703868674194879,0.071799982694269,0.278565453785017,0.418750000000002])
scatter(data_SLPIVASL_19(:,1), data_SLPIVASL_19(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_SLPIVASL_19', 100, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 0.5])
axis equal
text(0.018018018018018,0.932432432432432,0,'(f)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
xlabel('Local min u$^{\prime}/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
xlim([-6 6])
ylim([-6 6])

%% du(next),du(min,max) distribution

duwrtglobalmaxmin_HotWT7 = zeros(size(VF_HotWT7.u,1),size(VF_HotWT7.u,2)-1);
du_HotWT7 = abs(diff(VF_HotWT7uprimelong,1,2));
for z = 1:size(VF_HotWT7.z,1)
    [peaks, peak_locs] = findpeaks(VF_HotWT7uprimelong(z,:)); % Find local maxima
    [inverted_peaks, min_locs] = findpeaks(-VF_HotWT7uprimelong(z,:)); % Find local minima
    valleys = -inverted_peaks;
    
    for c = 1:size(VF_HotWT7.u,2)-1
        [~,i] = find(c <= peak_locs, 1, 'first');
        [~,j] = find(c <= min_locs, 1, 'first');

        if peak_locs(1,1) < min_locs(1,1)
            if c <= peak_locs(1,1)
                % Direct slope computation (signed)
                duwrtglobalmaxmin_HotWT7(z,c) = abs(VF_HotWT7uprimelong(z,peak_locs(1,i)) - VF_HotWT7uprimelong(z,c));
            elseif ~isempty(i) || ~isempty(j)
                % Compute both signed slopes
                slope_vals = [
                    VF_HotWT7uprimelong(z,peak_locs(1,i-1)) - VF_HotWT7uprimelong(z,c),...
                    VF_HotWT7uprimelong(z,min_locs(1,j)) - VF_HotWT7uprimelong(z,c)
                ];
                
                % Select signed slope corresponding to min(abs())
                [~, idx] = min(abs(slope_vals), [], 'includenan');
                duwrtglobalmaxmin_HotWT7(z,c) = abs(slope_vals(idx));
            else
                duwrtglobalmaxmin_HotWT7(z,c) = abs(VF_HotWT7uprimelong(z,min_locs(1,end)) - VF_HotWT7uprimelong(z,c));
            end
        else
            if c <= min_locs(1,1)
                duwrtglobalmaxmin_HotWT7(z,c) = abs(VF_HotWT7uprimelong(z,min_locs(1,j)) - VF_HotWT7uprimelong(z,c));
            elseif ~isempty(i) || ~isempty(j)
                slope_vals = [
                    VF_HotWT7uprimelong(z,peak_locs(1,i)) - VF_HotWT7uprimelong(z,c),...
                    VF_HotWT7uprimelong(z,min_locs(1,j-1)) - VF_HotWT7uprimelong(z,c)
                ];

                % Select signed slope corresponding to min(abs())
                [~, idx] = min(abs(slope_vals), [], 'includenan');
                duwrtglobalmaxmin_HotWT7(z,c) = abs(slope_vals(idx));
            else
                duwrtglobalmaxmin_HotWT7(z,c) = abs(VF_HotWT7uprimelong(z,peak_locs(1,end)) - VF_HotWT7uprimelong(z,c));
            end
        end
    end
end

figure
plot(duwrtglobalmaxmin_HotWT7(1,:),du_HotWT7(1,:),'k.')


duwrtglobalmaxmin_HotWT10 = zeros(size(VF_HotWT10.u,1),size(VF_HotWT10.u,2)-1);
du_HotWT10 = abs(diff(VF_HotWT10uprimelong,1,2));
for z = 1:size(VF_HotWT10.z,1)
    [peaks, peak_locs] = findpeaks(VF_HotWT10uprimelong(z,:)); % Find local maxima
    [inverted_peaks, min_locs] = findpeaks(-VF_HotWT10uprimelong(z,:)); % Find local minima
    valleys = -inverted_peaks;
    
    for c = 1:size(VF_HotWT7.u,2)-1
        [~,i] = find(c <= peak_locs, 1, 'first');
        [~,j] = find(c <= min_locs, 1, 'first');

        if peak_locs(1,1) < min_locs(1,1)
            if c <= peak_locs(1,1)
                % Direct slope computation (signed)
                duwrtglobalmaxmin_HotWT10(z,c) = abs(VF_HotWT10uprimelong(z,peak_locs(1,i)) - VF_HotWT10uprimelong(z,c));
            elseif ~isempty(i) || ~isempty(j)
                % Compute both signed slopes
                slope_vals = [
                    VF_HotWT10uprimelong(z,peak_locs(1,i-1)) - VF_HotWT10uprimelong(z,c),...
                    VF_HotWT10uprimelong(z,min_locs(1,j)) - VF_HotWT10uprimelong(z,c)
                ];
                
                % Select signed slope corresponding to min(abs())
                [~, idx] = min(abs(slope_vals), [], 'includenan');
                duwrtglobalmaxmin_HotWT10(z,c) = abs(slope_vals(idx));
            else
                duwrtglobalmaxmin_HotWT10(z,c) = abs(VF_HotWT10uprimelong(z,min_locs(1,end)) - VF_HotWT10uprimelong(z,c));
            end
        else
            if c <= min_locs(1,1)
                duwrtglobalmaxmin_HotWT10(z,c) = abs(VF_HotWT10uprimelong(z,min_locs(1,j)) - VF_HotWT10uprimelong(z,c));
            elseif ~isempty(i) || ~isempty(j)
                slope_vals = [
                    VF_HotWT10uprimelong(z,peak_locs(1,i)) - VF_HotWT10uprimelong(z,c),...
                    VF_HotWT10uprimelong(z,min_locs(1,j-1)) - VF_HotWT10uprimelong(z,c)
                ];

                % Select signed slope corresponding to min(abs())
                [~, idx] = min(abs(slope_vals), [], 'includenan');
                duwrtglobalmaxmin_HotWT10(z,c) = abs(slope_vals(idx));
            else
                duwrtglobalmaxmin_HotWT10(z,c) = abs(VF_HotWT10uprimelong(z,peak_locs(1,end)) - VF_HotWT10uprimelong(z,c));
            end
        end
    end
end

figure
plot(duwrtglobalmaxmin_HotWT10(10,:),du_HotWT10(10,:),'k.')


duwrtglobalmaxmin_SLPIVASL = zeros(size(VF_SLPIVASL.u,1),size(VF_SLPIVASL.u,3)-1);
du_SLPIVASL = abs(diff(VF_SLPIVASLuprimelong(:,25,:),1,3));
for z = 1:size(VF_SLPIVASL.z,1)
    [peaks, peak_locs] = findpeaks(reshape(VF_SLPIVASLuprimelong(z,25,:),1,[])); % Find local maxima
    [inverted_peaks, min_locs] = findpeaks(reshape(-VF_SLPIVASLuprimelong(z,25,:),1,[])); % Find local minima
    valleys = -inverted_peaks;
    
    for c = 1:size(VF_SLPIVASL.u,3)-1
        [~,i] = find(c <= peak_locs, 1, 'first');
        [~,j] = find(c <= min_locs, 1, 'first');

        if peak_locs(1,1) < min_locs(1,1)
            if c <= peak_locs(1,1)
                % Direct slope computation (signed)
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(VF_SLPIVASLuprimelong(z,25,peak_locs(1,i)) - VF_SLPIVASLuprimelong(z,25,c));
            elseif ~isempty(i) || ~isempty(j)
                % Compute both signed slopes
                slope_vals = [
                    VF_SLPIVASLuprimelong(z,25,peak_locs(1,i-1)) - VF_SLPIVASLuprimelong(z,25,c),...
                    VF_SLPIVASLuprimelong(z,25,min_locs(1,j)) - VF_SLPIVASLuprimelong(z,25,c)
                ];
                
                % Select signed slope corresponding to min(abs())
                [~, idx] = min(abs(slope_vals), [], 'includenan');
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(slope_vals(idx));
            elseif peak_locs(1,end) < min_locs(1,end)
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(VF_SLPIVASLuprimelong(z,25,min_locs(1,end)) - VF_SLPIVASLuprimelong(z,25,c));
            else
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(VF_SLPIVASLuprimelong(z,25,peak_locs(1,end)) - VF_SLPIVASLuprimelong(z,25,c));
            end
        else
            if c <= min_locs(1,1)
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(VF_SLPIVASLuprimelong(z,25,min_locs(1,j)) - VF_SLPIVASLuprimelong(z,25,c));
            elseif ~isempty(i) || ~isempty(j)
                slope_vals = [
                    VF_SLPIVASLuprimelong(z,25,peak_locs(1,i)) - VF_SLPIVASLuprimelong(z,25,c),...
                    VF_SLPIVASLuprimelong(z,25,min_locs(1,j-1)) - VF_SLPIVASLuprimelong(z,25,c)
                ];

                % Select signed slope corresponding to min(abs())
                [~, idx] = min(abs(slope_vals), [], 'includenan');
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(slope_vals(idx));
            elseif peak_locs(1,end) > min_locs(1,end)
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(VF_SLPIVASLuprimelong(z,25,peak_locs(1,end)) - VF_SLPIVASLuprimelong(z,25,c));
            else
                duwrtglobalmaxmin_SLPIVASL(z,c) = abs(VF_SLPIVASLuprimelong(z,25,min_locs(1,end)) - VF_SLPIVASLuprimelong(z,25,c));
            end
        end
    end
end

figure
plot(duwrtglobalmaxmin_SLPIVASL(10,:),du_SLPIVASL(10,:),'k.')

%% Saving & Loading

save('duwrtglobalmaxmin_HotWT7.mat','duwrtglobalmaxmin_HotWT7')
save('duwrtglobalmaxmin_HotWT10.mat','duwrtglobalmaxmin_HotWT10')
save('duwrtglobalmaxmin_SLPIVASL.mat','duwrtglobalmaxmin_SLPIVASL')

load('duwrtglobalmaxmin_HotWT7.mat')
load('duwrtglobalmaxmin_HotWT10.mat')
load('duwrtglobalmaxmin_SLPIVASL.mat')

%% Ploting

nbins = 1000;
data_HotWT7_1 = [duwrtglobalmaxmin_HotWT7(1,:)'/(VF_HotWT7.u_tau),...
    du_HotWT7(1,:)'/(VF_HotWT7.u_tau*res_WT7(1,1))]; % Example correlated non-Gaussian data

Y_min = min(data_HotWT7_1(:,2));
if Y_min <= 0
    Y_min = min(data_HotWT7_1(data_HotWT7_1(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_HotWT7_1(:,1)), max(data_HotWT7_1(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_HotWT7_1(:,2))), nbins+1);

% Compute 2D histogram
[counts, X_edges, Y_edges] = histcounts2(data_HotWT7_1(:,1), data_HotWT7_1(:,2), X_edges, Y_edges);

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (log-spaced in Y)
dx = mean(diff(X_edges)); % Bin width in X
dy = diff(log(Y_edges));  % Logarithmic bin width in Y
bin_area = dx * mean(dy); % Area of each bin (average dy for simplicity)

%  Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT7_1_normalized = counts / (sum(counts(:)) * bin_area);

% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_HotWT7_1_normalized(:)) * bin_area;
disp(['Renormalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

figure;

set(gcf,'Position',[653,546,806,584])
axes('Position',[0.08030112923463,0.54639175257732,0.265997490589711,0.36426116838488])
scatter(data_HotWT7_1(:,1), data_HotWT7_1(:,2), 5, [0.0,1,0.0],'filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT7_1_normalized', 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 1])
hcb2 = colorbar;
hcb2.Location = "north";
hcb2.Position = [0.097867001254705,0.924460431654677,0.887076537013795,0.025903057081893];
title(hcb2, 'p.d.f.', 'Interpreter', 'Latex', 'FontSize', 20, 'position',...
    [-32.62500000000094,-4.775,0])
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation = 'out';
axis square
% xlabel('Local min u$^{\prime}$','Interpreter','latex','FontSize',15)
ylabel('$\left|\mathrm{du}^{\prime}_{\mathrm{Neighbour}}\right|/(\Delta\mathrm{x}_{\mathrm{data}}\mathrm{u}_{\tau})$','Interpreter','latex','FontSize',15,...
    'position',[-1.034593624888726,3.162301964601093,-1])
text(0.83877273499915,0.097526772055074,0,'(a)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlim([0 6])
ylim([10^-3 10^4])

data_HotWT7_9 = [duwrtglobalmaxmin_HotWT7(9,:)'/(VF_HotWT7.u_tau),...
    du_HotWT7(9,:)'/(VF_HotWT7.u_tau*res_WT7(9,1))]; % Example correlated non-Gaussian data


Y_min = min(data_HotWT7_9(:,2));
if Y_min <= 0
    Y_min = min(data_HotWT7_9(data_HotWT7_9(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_HotWT7_9(:,1)), max(data_HotWT7_9(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_HotWT7_9(:,2))), nbins+1);

% Compute 2D histogram
[counts, X_edges, Y_edges] = histcounts2(data_HotWT7_9(:,1), data_HotWT7_9(:,2), X_edges, Y_edges);

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (log-spaced in Y)
dx = mean(diff(X_edges)); % Bin width in X
dy = diff(log(Y_edges));  % Logarithmic bin width in Y
bin_area = dx * mean(dy); % Area of each bin (average dy for simplicity)

%  Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT7_9_normalized = counts / (sum(counts(:)) * bin_area);


% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_HotWT7_9_normalized(:)) * bin_area;
disp(['Renormalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

axes('Position',[0.08030112923463,0.104810996563568,0.265997490589711,0.364261168384879])
scatter(data_HotWT7_9(:,1), data_HotWT7_9(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT7_9_normalized', 600, 'LineWidth', 3); % Contour of PDF
% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 1])
axis square
text(0.83877273499915,0.097526772055074,0,'(b)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
xlabel('$\left|\mathrm{du}^{\prime}_{\mathrm{min}\{\mathrm{Local\:max, Local\:min}\}}\right|/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
ylabel('$\left|\mathrm{du}^{\prime}_{\mathrm{Neighbour}}\right|/(\Delta\mathrm{x}_{\mathrm{data}}\mathrm{u}_{\tau})$','Interpreter','latex','FontSize',15,...
    'position',[-1.034593624888726,3.162301964601093,-1])
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlim([0 6])
ylim([10^-3 10^4])


data_HotWT10_2 = [duwrtglobalmaxmin_HotWT10(2,:)'/(VF_HotWT10.u_tau),...
    du_HotWT10(2,:)'/(VF_HotWT10.u_tau*res_WT10(2,1))]; % Example correlated non-Gaussian data

Y_min = min(data_HotWT10_2(:,2));
if Y_min <= 0
    Y_min = min(data_HotWT10_2(data_HotWT10_2(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_HotWT10_2(:,1)), max(data_HotWT10_2(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_HotWT10_2(:,2))), nbins+1);

% Compute 2D histogram
[counts, X_edges, Y_edges] = histcounts2(data_HotWT10_2(:,1), data_HotWT10_2(:,2), X_edges, Y_edges);

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (log-spaced in Y)
dx = mean(diff(X_edges)); % Bin width in X
dy = diff(log(Y_edges));  % Logarithmic bin width in Y
bin_area = dx * mean(dy); % Area of each bin (average dy for simplicity)

%  Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT10_2_normalized = counts / (sum(counts(:)) * bin_area);


% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_HotWT10_2_normalized(:)) * bin_area;
disp(['Renormalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

axes('Position',[0.401505646173145,0.54639175257732,0.265997490589709,0.36426116838488])
scatter(data_HotWT10_2(:,1), data_HotWT10_2(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT10_2_normalized', 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 1])
axis square
text(0.83877273499915,0.097526772055074,0,'(c)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlim([0 6])
ylim([10^-3 10^4])

data_HotWT10_9 = [duwrtglobalmaxmin_HotWT10(9,:)'/(VF_HotWT10.u_tau),...
    du_HotWT10(9,:)'/(VF_HotWT10.u_tau*res_WT10(9,1))]; % Example correlated non-Gaussian data

Y_min = min(data_HotWT10_9(:,2));
if Y_min <= 0
    Y_min = min(data_HotWT10_9(data_HotWT10_9(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_HotWT10_9(:,1)), max(data_HotWT10_9(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_HotWT10_9(:,2))), nbins+1);

% Compute 2D histogram
[counts, X_edges, Y_edges] = histcounts2(data_HotWT10_9(:,1), data_HotWT10_9(:,2), X_edges, Y_edges);

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (log-spaced in Y)
dx = mean(diff(X_edges)); % Bin width in X
dy = diff(log(Y_edges));  % Logarithmic bin width in Y
bin_area = dx * mean(dy); % Area of each bin (average dy for simplicity)

%  Normalize histogram to create probability density function (PDF)
pdf_hist_HotWT10_9_normalized = counts / (sum(counts(:)) * bin_area);

% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_HotWT10_9_normalized(:)) * bin_area;
disp(['Renormalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

axes('Position',[0.401505646173145,0.104810996563567,0.265997490589709,0.364261168384879])
scatter(data_HotWT10_9(:,1), data_HotWT10_9(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_HotWT10_9_normalized', 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 1])
axis square
text(0.83877273499915,0.097526772055074,0,'(d)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
xlabel('$\left|\mathrm{du}^{\prime}_{\mathrm{min}\{\mathrm{Local\:max, Local\:min}\}}\right|/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlim([0 6])
ylim([10^-3 10^4])

data_SLPIVASL_7 = [duwrtglobalmaxmin_SLPIVASL(7,:)'/(VF_SLPIVASL.u_tau),...
    du_SLPIVASL(7,:)'/(VF_SLPIVASL.u_tau*res_ASL(7,25))]; % Example correlated non-Gaussian data

Y_min = min(data_SLPIVASL_7(:,2));
if Y_min <= 0
    Y_min = min(data_SLPIVASL_7(data_SLPIVASL_7(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_SLPIVASL_7(:,1)), max(data_SLPIVASL_7(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_SLPIVASL_7(:,2))), nbins+1);

% Compute 2D histogram
[counts, X_edges, Y_edges] = histcounts2(data_SLPIVASL_7(:,1), data_SLPIVASL_7(:,2), X_edges, Y_edges);

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (log-spaced in Y)
dx = mean(diff(X_edges)); % Bin width in X
dy = diff(log(Y_edges));  % Logarithmic bin width in Y
bin_area = dx * mean(dy); % Area of each bin (average dy for simplicity)

%  Normalize histogram to create probability density function (PDF)
pdf_hist_SLPIVASL_7_log_normalized = counts / (sum(counts(:)) * bin_area);

% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_SLPIVASL_7_log_normalized(:)) * bin_area;
disp(['Renormalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

axes('Position',[0.718974068389199,0.54639175257732,0.265997490589708,0.36426116838488])
scatter(data_SLPIVASL_7(:,1), data_SLPIVASL_7(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_SLPIVASL_7_log_normalized', 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 1])
axis square
text(0.83877273499915,0.097526772055074,0,'(e)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlim([0 6])
ylim([10^-3 10^4])


data_SLPIVASL_19 = [duwrtglobalmaxmin_SLPIVASL(19,:)'/(VF_SLPIVASL.u_tau),...
    du_SLPIVASL(19,:)'/(VF_SLPIVASL.u_tau*res_ASL(19,25))]; % Example correlated non-Gaussian data

Y_min = min(data_SLPIVASL_19(:,2));
if Y_min <= 0
    Y_min = min(data_SLPIVASL_19(data_SLPIVASL_19(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_SLPIVASL_19(:,1)), max(data_SLPIVASL_19(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_SLPIVASL_19(:,2))), nbins+1);

% Compute 2D histogram
[counts, X_edges, Y_edges] = histcounts2(data_SLPIVASL_19(:,1), data_SLPIVASL_19(:,2), X_edges, Y_edges);

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Compute bin width (log-spaced in Y)
dx = mean(diff(X_edges)); % Bin width in X
dy = diff(log(Y_edges));  % Logarithmic bin width in Y
bin_area = dx * mean(dy); % Area of each bin (average dy for simplicity)

%  Normalize histogram to create probability density function (PDF)
pdf_hist_SLPIVASL_19_log_normalized = counts / (sum(counts(:)) * bin_area);

% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_SLPIVASL_19_log_normalized(:)) * bin_area;
disp(['Renormalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

axes('Position',[0.718974068389199,0.104810996563568,0.265997490589708,0.364261168384879])
scatter(data_SLPIVASL_19(:,1), data_SLPIVASL_19(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, pdf_hist_SLPIVASL_19_log_normalized', 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([0 1])
axis square
text(0.83877273499915,0.097526772055074,0,'(f)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
xlabel('$\left|\mathrm{du}^{\prime}_{\mathrm{min}\{\mathrm{Local\:max, Local\:min}\}}\right|/\mathrm{u}_{\tau}$','Interpreter','latex','FontSize',15)
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlim([0 6])
ylim([10^-3 10^4])

%%

figure
plot((0:0.0001:1-0.0001).*mean(VF_HotWT7.u(1,:),2),VF_HotWT7uprimelong(1,1:1:10000),'b.')

% figure
% plot((0:1/120:1-1/120).*mean(VF_SLPIVASL.u(1,25,:),3),reshape(VF_SLPIVASL.u(1,25,1:1:120),[],1))
% 
% 
% figure
% plot((0:1/120:1-1/120).*mean(VF_SLPIVASL.u(2,25,:),3),reshape(VF_SLPIVASL.u(2,25,1:1:120),[],1))

figure
for z = 1:size(VF_SLPIVASL.z,1)

plot((0:1/120:1-1/120).*mean(VF_SLPIVASL.u(z,25,:),3),reshape(VF_SLPIVASL.u(z,25,1:1:120),[],1))
hold on
    if z>=25
        [r,lags] = xcorr(reshape(VF_SLPIVASL.uprime(z-24,25,:),[],1),reshape(VF_SLPIVASL.uprime(z,25,:),[],1),'coeff');
        zero_lag_index = find(lags == 0); % Find index where lag is zero
        zero_lag_corr = r(zero_lag_index)
    end
end


for z = 1:size(VF_HotWT7.z,1)
    if z~=1
        [r,lags] = xcorr(VF_HotWT7.uprime(z-1,:),VF_HotWT7.uprime(z,:),'coeff');
        zero_lag_index = find(lags == 0); % Find index where lag is zero
        zero_lag_corr = r(zero_lag_index)
    end
end

%%
% Define parameters
num_rows = 50;   % Number of rows (iterations)
num_cols = 1000; % Length of each signal
alpha = 0.99;    % Desired correlation with the previous row

% Initialize matrix
X = zeros(num_rows, num_cols);

% Generate the first row with random values
X(1, :) = randn(1, num_cols); % First random signal

% Generate subsequent rows with correlation of 0.99 to the previous row
for i = 2:num_rows
    noise = randn(1, num_cols); % Independent random noise
    X(i, :) = alpha * X(i-1, :) + sqrt(1 - alpha^2) * noise; % Correlated signal
end

% Verify correlation between rows
correlation_matrix = corrcoef(X'); % Compute row-wise correlation matrix
disp('Correlation matrix (subset of first 10 rows):');
disp(correlation_matrix(1:10, 1:10)); % Display a subset for verification

% Plot a few rows for visualization
figure;
hold on;
for i = 1:5  % Plot first 5 rows
    plot(X(i, :), 'DisplayName', ['Row ', num2str(i)]);
end
legend show;
title('First 5 Rows of Generated Correlated Signals');
xlabel('Time');
ylabel('Amplitude');
grid on;
