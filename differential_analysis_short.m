Data3 = load('VF_HotWT7.mat');
VF_HotWT7 = Data3.VF_HotWT7;
clear('Data3')


Data6 = load('VF_HotWT10.mat');
VF_HotWT10 = Data6.VF_HotWT10;
clear('Data6')


Data8 = load('VF_SLPIVASL.mat');
VF_SLPIVASL = Data8.VF_SLPIVASL;
clear('Data8')

Data9 = load('VF_SonicASL.mat');
VF_SonicASL = Data9.VF_SonicASL;
clear('Data9')
%% Long decomposition
VF_HotWT7uprimelong = VF_HotWT7.u - mean(VF_HotWT7.u,2);
VF_HotWT7wprimelong = VF_HotWT7.w;
VF_HotWT10uprimelong = VF_HotWT10.u - mean(VF_HotWT10.u,2);
VF_HotWT10wprimelong = VF_HotWT10.w;
VF_SLPIVASLuprimelong = VF_SLPIVASL.uprime;
VF_SLPIVASLwprimelong = VF_SLPIVASL.w;

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
res_WT7 = X_domain_WT7(:,2)-X_domain_WT7(:,1);

T_domain_WT10 = (0:1:size(VF_HotWT10.u,2)-1)*1/VF_HotWT10.fs;
X_domain_WT10 = mean(VF_HotWT10.u,2)*T_domain_WT10;
res_WT10 = X_domain_WT10(:,2)-X_domain_WT10(:,1);

T_domain_ASL = (0:1:size(VF_SLPIVASL.u,3)-1)*1/VF_SLPIVASL.fs;
T_domain_ASL = reshape(T_domain_ASL,1,1,[]);
X_domain_ASL = mean(VF_SLPIVASL.u,3).*T_domain_ASL;
res_ASL = X_domain_ASL(:,:,2)-X_domain_ASL(:,:,1);


VF_HotWT7uprimeshort = cell(size(VF_HotWT7.z,1), 1);
VF_HotWT7wprimeshort = cell(size(VF_HotWT7.z,1), 1);
VF_HotWT10uprimeshort = cell(size(VF_HotWT10.z,1), 1);
VF_HotWT10wprimeshort = cell(size(VF_HotWT10.z,1), 1);
VF_SLPIVASLuprimeshort = cell(size(VF_SLPIVASL.z,1), 1);
VF_SLPIVASLwprimeshort = cell(size(VF_SLPIVASL.z,1), 1);

for i = 1:size(VF_SLPIVASLuprimeshort,1)
    VF_SLPIVASLuprimeshort{i} = cell(1, size(VF_SLPIVASL.u,2));
    VF_SLPIVASLwprimeshort{i} = cell(1, size(VF_SLPIVASL.u,2));
end

C = zeros();
% figure
% set(gcf,'Position',[814,534,806,379])
% axes('Position',[0.064516129032258,0.12664907651715,0.915632754342432,0.798350923482849])

l0WT7 = 0.7;%0.3
l0WT10 = 0.7;%0.3
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
        VF_HotWT7uprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2)+1e-6*rand(1,unitindex_WT7);%%added random value useful when I delete zero arrays
        
        VF_HotWT7wprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7) = ...
            VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)-...
            mean(VF_HotWT7.w(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),2)+1e-6*rand(1,unitindex_WT7);
        
        
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/l0WT7,...
%             VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
%         hold on
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/l0WT7,...
%             VF_HotWT7uprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
% 
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/l0WT7,...
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

%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/(c0WT7*lambda_T_WT7(z,1)),...
%         VF_HotWT7.u(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
%         hold on
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/(c0WT7*lambda_T_WT7(z,1)),...
%             VF_HotWT7uprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7),'linewidth',2)
% 
%         plot(X_domain_WT7(z,(s-1)*unitindex_WT7+1:(s)*unitindex_WT7)/(c0WT7*lambda_T_WT7(z,1)),...
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
        
        R = corrcoef(VF_HotWT7uprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7)',...
            VF_HotWT7wprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7)');
        C(z,s) = R(1,2);
    end
    %Add ks density
    
end


for z=1:size(VF_HotWT7.z,1)
    
    VF_HotWT7.z(z)/VF_HotWT7.delta
    R = corrcoef(VF_HotWT7uprimeshort{z},VF_HotWT7wprimeshort{z})
    
end

for z = 1: size(VF_HotWT10.z,1)
    
%     ratio_WT10 = X_domain_WT10(z,:) / (c0WT10*lambda_T_WT10(z,1));
    ratio_WT10 = X_domain_WT10(z,:) / l0WT10;
    unitindex_WT10 = find(ratio_WT10 >= 1, 1, 'first');
    unitindex_WT10 = unitindex_WT10 - 1;
    
    for s = 1 : floor(size(VF_HotWT10.u(z,:),2)/unitindex_WT10)
        VF_HotWT10uprimeshort{z}((s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.u(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2)+1e-6*rand(1,unitindex_WT10);
        
        VF_HotWT10wprimeshort{z}((s-1)*unitindex_WT10+1:(s)*unitindex_WT10) = ...
            VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10)-...
            mean(VF_HotWT10.w(z,(s-1)*unitindex_WT10+1:(s)*unitindex_WT10),2)+1e-6*rand(1,unitindex_WT10);
    end
    
end



for z=1:size(VF_HotWT10.z,1)
    
    VF_HotWT10.z(z)/VF_HotWT10.delta
    R = corrcoef(VF_HotWT10uprimeshort{z},VF_HotWT10wprimeshort{z})
    
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
            VF_SLPIVASLuprimeshort{z}{j}((s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)) = ...
                VF_SLPIVASL.u(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j))-...
                mean(VF_SLPIVASL.u(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)),3)+1e-6*rand(1,1,unitindex_ASL(j));
            
            VF_SLPIVASLwprimeshort{z}{j}((s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)) = ...
                VF_SLPIVASL.w(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j))-...
                mean(VF_SLPIVASL.w(z,j,(s-1)*unitindex_ASL(1,j)+1:(s)*unitindex_ASL(1,j)),3)+1e-6*rand(1,1,unitindex_ASL(j));
            
        end
    end
end



for z=1:size(VF_SLPIVASL.z,1)
    
    VF_SLPIVASL.z(z)/VF_SLPIVASL.delta
    R = corrcoef(VF_SLPIVASLuprimeshort{z}{25},VF_SLPIVASLwprimeshort{z}{25})
    
end
%% min-max distribution short

value_pairs_cell_HotWT7 = cell(size(VF_HotWT7.z,1), 1);

for z = 1:size(VF_HotWT7.z,1)
    
%     ratio_WT7= X_domain_WT7(z,:) / (c0WT7*lambda_T_WT7(z,1));
    ratio_WT7= X_domain_WT7(z,:) / l0WT7;
    unitindex_WT7 = find(ratio_WT7 >= 1, 1, 'first');
    unitindex_WT7 = unitindex_WT7 - 1;
    % Initialize an empty matrix for storing value pairs at this elevation
    value_pairs = [];
    
    
    for s = 1 : floor(size(VF_HotWT7uprimeshort{z},2)/unitindex_WT7)
        % Find local maxima and minima for the current elevation
        [peaks, peak_locs] = findpeaks(VF_HotWT7uprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7)); % Local maxima from long should be changed to short
        [inverted_peaks, min_locs] = findpeaks(-VF_HotWT7uprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7)); % Local minima
        valleys = -inverted_peaks;
        
        % Check if there are any maxima or minima found
        if isempty(peak_locs) || isempty(min_locs)
%             disp(['No local maxima or minima found at elevation z = ', num2str(z)]);
            continue
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
    end
    % Store the value pairs in the cell array for this elevation
    value_pairs_cell_HotWT7{z} = value_pairs;
    
    
end


value_pairs_cell_HotWT10 = cell(size(VF_HotWT10.z,1), 1);

for z = 1:size(VF_HotWT10.z,1)
    
%     ratio_WT10 = X_domain_WT10(z,:) / (c0WT10*lambda_T_WT10(z,1));
    ratio_WT10 = X_domain_WT10(z,:) / l0WT10;
    unitindex_WT10 = find(ratio_WT10 >= 1, 1, 'first');
    unitindex_WT10 = unitindex_WT10 - 1;
    % Initialize an empty matrix for storing value pairs at this elevation
    value_pairs = [];
    
    for s = 1 : floor(size(VF_HotWT10uprimeshort{z},2)/unitindex_WT10)
        % Find local maxima and minima for the current elevation
        [peaks, peak_locs] = findpeaks(VF_HotWT10uprimeshort{z}((s-1)*unitindex_WT10+1:(s)*unitindex_WT10)); % Local maxima
        [inverted_peaks, min_locs] = findpeaks(-VF_HotWT10uprimeshort{z}((s-1)*unitindex_WT10+1:(s)*unitindex_WT10)); % Local minima
        valleys = -inverted_peaks;
        
        
        % Check if there are any maxima or minima found
        if isempty(peak_locs) || isempty(min_locs)
%             disp(['No local maxima or minima found at elevation z = ', num2str(z)]);
            continue
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
    end
    
    
    % Store the value pairs in the cell array for this elevation
    value_pairs_cell_HotWT10{z} = value_pairs;
    
    
end



value_pairs_cell_SLPIVASL = cell(size(VF_SLPIVASL.z,1), 1);

for z = 1:size(VF_SLPIVASL.z,1)
    
    
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
    % Initialize an empty matrix for storing value pairs at this elevation
    value_pairs = [];
    
    for k = 25%1:size(unitindex_ASL,2)
        
        for s = 1 : floor(size(VF_SLPIVASLuprimeshort{z}{k},2)/unitindex_ASL(1,k))
            % Find local maxima and minima for the current elevation
            [peaks, peak_locs] = findpeaks(reshape(VF_SLPIVASLuprimeshort{z}{k}((s-1)*unitindex_ASL(1,k)+1:(s)*unitindex_ASL(1,k)),1,[])); % Local maxima
            [inverted_peaks, min_locs] = findpeaks(reshape(-VF_SLPIVASLuprimeshort{z}{k}((s-1)*unitindex_ASL(1,k)+1:(s)*unitindex_ASL(1,k)),1,[]));
            valleys = -inverted_peaks;

            % Check if there are any maxima or minima found
            if isempty(peak_locs) || isempty(min_locs)
%                 disp(['No local maxima or minima found at elevation z = ', num2str(z)]);
                continue
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
        end
    end
    
    % Store the value pairs in the cell array for this elevation
    value_pairs_cell_SLPIVASL{z} = value_pairs;
    
    
end


%% Saving & Loading

save('value_pairs_cell_SLPIVASL.mat','value_pairs_cell_SLPIVASL')


load('value_pairs_cell_SLPIVASL.mat')


%% ploting ( WT7,WT10)
    
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

set(gcf,'Position',[904,437,755,487])
axes('Position',[0.075496688741721,0.077499427040178,0.415894039735099,0.725173210161663])
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
hcb2.Position = [0.109181141439206,0.901022931654678,0.871998281396428,0.033366208616816];
title(hcb2, 'j.p.d.f.', 'Interpreter', 'Latex', 'FontSize', 20, 'position',...
    [-29.625000000000938,-4.775,0])
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation = 'out';
axis equal

legend(sprintf('Pairwise local min-max $\\mathrm{u}^{\\prime}$ H-W(m1) z$/\\delta$ = %.2f',VF_HotWT7.z(1)/VF_HotWT7.delta),...
    'Histogram density estimation',...
    'Interpreter','latex','FontSize',10,'Position',...
    [0.066726049082404,0.803425808916316,0.428792657254608,0.074948663339477],...
    'Numcolumns',1,'Orientation','vertical','color','none');
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
xlabel('Local min sub--window u$^{\prime}/u_{\tau}$','Interpreter','latex','FontSize',16.5)
ylabel('Local max sub--window u$^{\prime}/u_{\tau}$','Interpreter','latex','FontSize',16.5)
text(0.018018018018018,0.932432432432432,0,'(a)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);

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


axes('Position',[0.577483443708607,0.077499427040178,0.415894039735101,0.725173210161663])
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
text(0.018018018018018,0.932432432432432,0,'(b)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',15);
legend(sprintf('Pairwise local min-max $\\mathrm{u}^{\\prime}$ H-W(m2) z$/\\delta$ = %.2f',VF_HotWT10.z(2)/VF_HotWT7.delta),...
    'Histogram density estimation',...
    'Interpreter','latex','FontSize',10,'Position',...
    [0.556963981211434,0.803425808916317,0.436396263195222,0.074948663339477],...
    'Numcolumns',1,'Orientation','vertical','color','none');
xlabel('Local min sub--window u$^{\prime}/u_{\tau}$','Interpreter','latex','FontSize',15)
ylabel('Local max sub--window u$^{\prime}/u_{\tau}$','Interpreter','latex','FontSize',15)
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
xlim([-6 6])
ylim([-6 6])



%% ploting (ALL: WT7,WT10,ASL)

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

%% du(next neighbor),du(min,max) distribution



duwrtlocalmaxmin_HotWT7 = cell(size(VF_HotWT7.z,1), 1);
duNN_HotWT7 = cell(size(VF_HotWT7.z,1), 1);

for z = 1:size(VF_HotWT7.z,1)

%     ratio_WT7= X_domain_WT7(z,:) / (c0WT7*lambda_T_WT7(z,1));
    ratio_WT7= X_domain_WT7(z,:) / l0WT7;
    unitindex_WT7 = find(ratio_WT7 >= 1, 1, 'first');
    unitindex_WT7 = unitindex_WT7 - 1;
    
    
    for s = 1 : floor(size(VF_HotWT7uprimeshort{z},2)/unitindex_WT7)
        u_temp = VF_HotWT7uprimeshort{z}((s-1)*unitindex_WT7+1:(s)*unitindex_WT7);
        [peaks, peak_locs] = findpeaks(u_temp); % Find local maxima
        [inverted_peaks, min_locs] = findpeaks(-u_temp); % Find local minima
        valleys = -inverted_peaks;
        if isempty(peak_locs) || isempty(min_locs)
%             disp(['No local maxima or minima found at elevation z = ', num2str(z)]);
            continue
        else
%             duNN_HotWT7{z}((s-1)*(unitindex_WT7-1)+1:(s)*(unitindex_WT7-1)) = abs(diff(u_temp));
            duNN_HotWT7{z} = [duNN_HotWT7{z},abs(diff(u_temp))];

            for c = 1:unitindex_WT7-1
                [~,i] = find(c <= peak_locs, 1, 'first');
                [~,j] = find(c <= min_locs, 1, 'first');
                
                if peak_locs(1,1) < min_locs(1,1)
                    if c <= peak_locs(1,1)
                        % Direct slope computation (signed)
%                         duwrtlocalmaxmin_HotWT7{z}((s-1)*(unitindex_WT7-1)+c) = abs(u_temp(1,peak_locs(1,i)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT7{z} = [duwrtlocalmaxmin_HotWT7{z},abs(u_temp(1,peak_locs(1,i)) - u_temp(1,c))];
                    elseif ~isempty(i) || ~isempty(j)
                        % Compute both signed slopes
                        slope_vals = [
                            u_temp(1,peak_locs(1,i-1)) - u_temp(1,c),...
                            u_temp(1,min_locs(1,j)) - u_temp(1,c)
                            ];
                        
                        % Select signed slope corresponding to min(abs())
                        [~, idx] = min(abs(slope_vals), [], 'includenan');
%                         duwrtlocalmaxmin_HotWT7{z}((s-1)*(unitindex_WT7-1)+c)  = abs(slope_vals(idx));
                        duwrtlocalmaxmin_HotWT7{z} = [duwrtlocalmaxmin_HotWT7{z},abs(slope_vals(idx))];
                    else
%                         duwrtlocalmaxmin_HotWT7{z}((s-1)*(unitindex_WT7-1)+c)  = abs(u_temp(1,min_locs(1,end)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT7{z} = [duwrtlocalmaxmin_HotWT7{z},abs(u_temp(1,min_locs(1,end)) - u_temp(1,c))];
                    end
                else
                    if c <= min_locs(1,1)
%                         duwrtlocalmaxmin_HotWT7{z}((s-1)*(unitindex_WT7-1)+c)  = abs(u_temp(1,min_locs(1,j)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT7{z} = [duwrtlocalmaxmin_HotWT7{z},abs(u_temp(1,min_locs(1,j)) - u_temp(1,c))];
                    elseif ~isempty(i) || ~isempty(j)
                        slope_vals = [
                            u_temp(1,peak_locs(1,i)) - u_temp(1,c),...
                            u_temp(1,min_locs(1,j-1)) - u_temp(1,c)
                            ];
                        
                        % Select signed slope corresponding to min(abs())
                        [~, idx] = min(abs(slope_vals), [], 'includenan');
%                         duwrtlocalmaxmin_HotWT7{z}((s-1)*(unitindex_WT7-1)+c)  = abs(slope_vals(idx));
                        duwrtlocalmaxmin_HotWT7{z}  = [duwrtlocalmaxmin_HotWT7{z},abs(slope_vals(idx))];
                    else
%                         duwrtlocalmaxmin_HotWT7{z}((s-1)*(unitindex_WT7-1)+c)  = abs(u_temp(1,peak_locs(1,end)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT7{z} = [duwrtlocalmaxmin_HotWT7{z},abs(u_temp(1,peak_locs(1,end)) - u_temp(1,c))];
                    end
                end
            end
        end
    end
end


figure
plot(duwrtlocalmaxmin_HotWT7{1},duNN_HotWT7{1},'k.')



duwrtlocalmaxmin_HotWT10 = cell(size(VF_HotWT10.z,1), 1);
duNN_HotWT10 = cell(size(VF_HotWT10.z,1), 1);

for z = 1:size(VF_HotWT10.z,1)

%     ratio_WT10 = X_domain_WT10(z,:) / (c0WT10*lambda_T_WT10(z,1));
    ratio_WT10= X_domain_WT10(z,:) / l0WT10;
    unitindex_WT10 = find(ratio_WT10 >= 1, 1, 'first');
    unitindex_WT10 = unitindex_WT10 - 1;
    
    
    for s = 1 : floor(size(VF_HotWT10uprimeshort{z},2)/unitindex_WT10)
        u_temp = VF_HotWT10uprimeshort{z}((s-1)*unitindex_WT10+1:(s)*unitindex_WT10);
        [peaks, peak_locs] = findpeaks(u_temp); % Find local maxima
        [inverted_peaks, min_locs] = findpeaks(-u_temp); % Find local minima
        valleys = -inverted_peaks;
        if isempty(peak_locs) || isempty(min_locs)
%             disp(['No local maxima or minima found at elevation z = ', num2str(z), ', index s = ', num2str(s)]);
            continue
        else
%             duNN_HotWT10{z}((s-1)*(unitindex_WT10-1)+1:(s)*(unitindex_WT10-1)) = abs(diff(u_temp));
            duNN_HotWT10{z}= [duNN_HotWT10{z},abs(diff(u_temp))];
            for c = 1:unitindex_WT10-1
                [~,i] = find(c <= peak_locs, 1, 'first');
                [~,j] = find(c <= min_locs, 1, 'first');
                
                if peak_locs(1,1) < min_locs(1,1)
                    if c <= peak_locs(1,1)
                        % Direct slope computation (signed)
%                         duwrtlocalmaxmin_HotWT10{z}((s-1)*(unitindex_WT10-1)+c) = abs(u_temp(1,peak_locs(1,i)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT10{z} = [duwrtlocalmaxmin_HotWT10{z},abs(u_temp(1,peak_locs(1,i)) - u_temp(1,c))];
                    elseif ~isempty(i) || ~isempty(j)
                        % Compute both signed slopes
                        slope_vals = [
                            u_temp(1,peak_locs(1,i-1)) - u_temp(1,c),...
                            u_temp(1,min_locs(1,j)) - u_temp(1,c)
                            ];
                        
                        % Select signed slope corresponding to min(abs())
                        [~, idx] = min(abs(slope_vals), [], 'includenan');
%                         duwrtlocalmaxmin_HotWT10{z}((s-1)*(unitindex_WT10-1)+c)  = abs(slope_vals(idx));
                        duwrtlocalmaxmin_HotWT10{z} = [duwrtlocalmaxmin_HotWT10{z},abs(slope_vals(idx))];
                    else
%                         duwrtlocalmaxmin_HotWT10{z}((s-1)*(unitindex_WT10-1)+c)  = abs(u_temp(1,min_locs(1,end)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT10{z} = [duwrtlocalmaxmin_HotWT10{z},abs(u_temp(1,min_locs(1,end)) - u_temp(1,c))];
                    end
                else
                    if c <= min_locs(1,1)
%                         duwrtlocalmaxmin_HotWT10{z}((s-1)*(unitindex_WT10-1)+c)  = abs(u_temp(1,min_locs(1,j)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT10{z} = [duwrtlocalmaxmin_HotWT10{z},abs(u_temp(1,min_locs(1,j)) - u_temp(1,c))];
                    elseif ~isempty(i) || ~isempty(j)
                        slope_vals = [
                            u_temp(1,peak_locs(1,i)) - u_temp(1,c),...
                            u_temp(1,min_locs(1,j-1)) - u_temp(1,c)
                            ];
                        
                        % Select signed slope corresponding to min(abs())
                        [~, idx] = min(abs(slope_vals), [], 'includenan');
%                         duwrtlocalmaxmin_HotWT10{z}((s-1)*(unitindex_WT10-1)+c)  = abs(slope_vals(idx));
                        duwrtlocalmaxmin_HotWT10{z} = [duwrtlocalmaxmin_HotWT10{z},abs(slope_vals(idx))];
                    else
%                         duwrtlocalmaxmin_HotWT10{z}((s-1)*(unitindex_WT10-1)+c)  = abs(u_temp(1,peak_locs(1,end)) - u_temp(1,c));
                        duwrtlocalmaxmin_HotWT10{z} = [duwrtlocalmaxmin_HotWT10{z},abs(u_temp(1,peak_locs(1,end)) - u_temp(1,c))];
                    end
                end
            end
        end
    end
end


figure
plot(duwrtlocalmaxmin_HotWT10{10},duNN_HotWT10{10},'k.')



duwrtlocalmaxmin_SLPIVASL = cell(size(VF_SLPIVASL.z,1), 1);
duNN_SLPIVASL = cell(size(VF_SLPIVASL.z,1), 1);

for i = 1:size(VF_SLPIVASLuprimeshort,1)
    duwrtlocalmaxmin_SLPIVASL{i} = cell(1, size(VF_SLPIVASL.u,2));
    duNN_SLPIVASL{i} = cell(1, size(VF_SLPIVASL.u,2));
end

for z = 1:size(VF_SLPIVASL.z,1)
    
    
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
    
    
    for k = 25%1:size(unitindex_ASL,2)
        counter = 1;
        for s = 1 : floor(size(VF_SLPIVASLuprimeshort{z}{k},2)/unitindex_ASL(1,k))
            % Find local maxima and minima for the current elevation
            u_temp = VF_SLPIVASLuprimeshort{z}{k}((s-1)*unitindex_ASL(1,k)+1:(s)*unitindex_ASL(1,k));
            [peaks, peak_locs] = findpeaks(u_temp); % Local maxima
            [inverted_peaks, min_locs] = findpeaks(-u_temp);
            valleys = -inverted_peaks;
            if isempty(peak_locs) || isempty(min_locs)
                disp(['No local maxima or minima found at elevation z = ', num2str(z), ', index s = ', num2str(s)]);
            else
                
                for c = 1:unitindex_ASL(1,k)-1
                    duNN_SLPIVASL{z}{k}(counter) = abs(u_temp(c+1)-u_temp(c));
                    [~,i] = find(c <= peak_locs, 1, 'first');
                    [~,j] = find(c <= min_locs, 1, 'first');
                    
                    if peak_locs(1,1) < min_locs(1,1)
                        if c <= peak_locs(1,1)
                            % Direct slope computation (signed)
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(u_temp(1,peak_locs(1,i)) - u_temp(1,c));
                            counter = counter+1;
                        elseif ~isempty(i) || ~isempty(j)
                            % Compute both signed slopes
                            slope_vals = [
                                u_temp(1,peak_locs(1,i-1)) - u_temp(1,c),...
                                u_temp(1,min_locs(1,j)) - u_temp(1,c)
                                ];
                            
                            % Select signed slope corresponding to min(abs())
                            [~, idx] = min(abs(slope_vals), [], 'includenan');
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(slope_vals(idx));
                            counter = counter+1;
                        elseif peak_locs(1,end) < min_locs(1,end)
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(u_temp(1,min_locs(1,end)) - u_temp(1,c));
                            counter = counter+1;
                        else
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(u_temp(1,peak_locs(1,end)) - u_temp(1,c));
                            counter = counter+1;
                        end
                    else
                        if c <= min_locs(1,1)
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(u_temp(1,min_locs(1,j)) - u_temp(1,c));
                            counter = counter+1;
                        elseif ~isempty(i) || ~isempty(j)
                            slope_vals = [
                                u_temp(1,peak_locs(1,i)) - u_temp(1,c),...
                                u_temp(1,min_locs(1,j-1)) - u_temp(1,c)
                                ];
                            
                            % Select signed slope corresponding to min(abs())
                            [~, idx] = min(abs(slope_vals), [], 'includenan');
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(slope_vals(idx));
                            counter = counter+1;
                        elseif peak_locs(1,end) > min_locs(1,end)
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(u_temp(1,peak_locs(1,end)) - u_temp(1,c));
                            counter = counter+1;
                        else
                            duwrtlocalmaxmin_SLPIVASL{z}{k}(counter)  = abs(u_temp(1,min_locs(1,end)) - u_temp(1,c));
                            counter = counter+1;
                        end
                    end
                    
                end
            end
        end
    end
    
end


figure
plot(duwrtlocalmaxmin_SLPIVASL{10}{25},duNN_SLPIVASL{10}{25},'k.')

%% Saving & Loading

VF_HotWT7 = struct(VF_HotWT7);
VF_HotWT10 = struct(VF_HotWT10);

save('SGFV_param_HotWT7_Long.mat','duwrtlocalmaxmin_HotWT7',...
    'value_pairs_cell_HotWT7','duNN_HotWT7','VF_HotWT7','T_domain_WT7',...
    'X_domain_WT7','res_WT7')


save('SGFV_param_HotWT10_Long.mat','duwrtlocalmaxmin_HotWT10',...
    'value_pairs_cell_HotWT10','duNN_HotWT10','VF_HotWT10','T_domain_WT10',...
    'X_domain_WT10','res_WT10')

save('duwrtlocalmaxmin_SLPIVASL.mat','duwrtlocalmaxmin_SLPIVASL')
save('duNN_SLPIVASL.mat','duNN_SLPIVASL')


load('SGFV_param_HotWT7.mat')
load('SGFV_param_HotWT10.mat')
load('duwrtlocalmaxmin_SLPIVASL.mat')




%% Ploting WT7,WT10

nbins = 1000;

data_HotWT7_1 = [duwrtlocalmaxmin_HotWT7{1}'/(VF_HotWT7.u_tau),...
    duNN_HotWT7{1}'/(VF_HotWT7.u_tau*res_WT7(1,1))]; % Example correlated non-Gaussian data

Y_min = min(data_HotWT7_1(:,2));
if Y_min <= 0
    Y_min = min(data_HotWT7_1(data_HotWT7_1(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_HotWT7_1(:,1)), max(data_HotWT7_1(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_HotWT7_1(:,2))), nbins+1);

% Compute 2D histogram
[pdf_hist_HotWT7_1_normalized, X_edges, Y_edges] = histcounts2(data_HotWT7_1(:,1), data_HotWT7_1(:,2), X_edges, Y_edges,'Normalization','pdf');

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);





% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_HotWT7_1_normalized .*(diff(X_edges)'* diff(Y_edges)),'all');
disp(['Normalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

figure;

set(gcf,'Position',[742,520,755,487])
axes('Position',[0.080794701986754,0.077499427040178,0.415894039735099,0.725173210161663])
scatter(data_HotWT7_1(:,1), data_HotWT7_1(:,2), 5, [0.0,1,0.0],'filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, log(pdf_hist_HotWT7_1_normalized'), 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([-10 1])
hcb2 = colorbar;
hcb2.Location = "north";
hcb2.Position = [0.184105960264901,0.901022931654678,0.805020482438283,0.033366208616816];
title(hcb2, 'log(j.p.d.f.)', 'Interpreter', 'Latex', 'FontSize', 20, 'position',...
    [-52.12500000000094,-4.775,0])
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation = 'out';
axis square
xlabel('$\min\{\left|\mathrm{du}^{\prime}_{\mathrm{NLmin}}\right|, \left|\mathrm{du}^{\prime}_{\mathrm{NLmax}}\right|\}/u_{\tau}$','Interpreter','latex','FontSize',15)
ylabel('$\left|\mathrm{du}^{\prime}_{\mathrm{NN}}\right|/(\Delta\mathrm{x}_{\mathrm{data}}u_{\tau})$','Interpreter','latex','FontSize',15,...
    'position',[-0.633319749858728,3.162301964601093,-1])
text(0.890629482986171,0.059820967464279,0,'(a)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',16.5);


legend(sprintf('$\\mathrm{Pairwise}\\ \\frac{\\min\\left\\{\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NLmin}}\\right|,\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NLmax}}\\right|\\right\\}}{u_{\\tau}},\\frac{\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NN}}\\right|}{\\Delta\\mathrm{x}_{\\mathrm{data}}u_{\\tau}}$ H-W(m1) $\\mathrm{z}/\\delta=%.2f$', VF_HotWT7.z(1)/VF_HotWT7.delta),...
    'Histogram density estimation',...
    'Interpreter','latex','FontSize',9.5,'Position',...
    [0.005868118044307,0.79777348666809,0.492230373635437,0.086253307835928],...
    'Numcolumns',1,'Orientation','vertical','color','none');
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlim([0 6])
ylim([10^-3 10^4])


data_HotWT10_2 = [duwrtlocalmaxmin_HotWT10{2}'/(VF_HotWT10.u_tau),...
    duNN_HotWT10{2}'/(VF_HotWT10.u_tau*res_WT10(2,1))]; % Example correlated non-Gaussian data

Y_min = min(data_HotWT10_2(:,2));
if Y_min <= 0
    Y_min = min(data_HotWT10_2(data_HotWT10_2(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_HotWT10_2(:,1)), max(data_HotWT10_2(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_HotWT10_2(:,2))), nbins+1);

% Compute 2D histogram
[pdf_hist_HotWT10_2_normalized, X_edges, Y_edges] = histcounts2(data_HotWT10_2(:,1), data_HotWT10_2(:,2), X_edges, Y_edges,'Normalization','pdf');

%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);



% Compute new total probability (should be ≈1 after renormalization)
P_final = sum(pdf_hist_HotWT10_2_normalized .*(diff(X_edges)'* diff(Y_edges)),'all');
disp(['Normalized Total Probability (should be ≈ 1): ', num2str(P_final)]);

axes('Position',[0.577483443708607,0.077499427040178,0.415894039735101,0.725173210161663])
scatter(data_HotWT10_2(:,1), data_HotWT10_2(:,2), 5, 'g','filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, log(pdf_hist_HotWT10_2_normalized'), 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([-10 1])
axis square
text(0.883358722260296,0.05931021154552,0,'(b)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',16.5);
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
legend(sprintf('$\\mathrm{Pairwise}\\ \\frac{\\min\\left\\{\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NLmin}}\\right|,\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NLmax}}\\right|\\right\\}}{u_{\\tau}},\\frac{\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NN}}\\right|}{\\Delta\\mathrm{x}_{\\mathrm{data}}u_{\\tau}}$ H-W(m2) $\\mathrm{z}/\\delta=%.2f$',VF_HotWT10.z(2)/VF_HotWT7.delta),...
    'Histogram density estimation',...
    'Interpreter','latex','FontSize',9.5,'Position',...
    [0.505205866388678,0.797773486668092,0.492230373635437,0.086253307835928],...
    'Numcolumns',1,'Orientation','vertical','color','none');
xlabel('$\min\{\left|\mathrm{du}^{\prime}_{\mathrm{NLmin}}\right|, \left|\mathrm{du}^{\prime}_{\mathrm{NLmax}}\right|\}/u_{\tau}$','Interpreter','latex','FontSize',15)
ylabel('$\left|\mathrm{du}^{\prime}_{\mathrm{NN}}\right|/(\Delta\mathrm{x}_{\mathrm{data}}u_{\tau})$','Interpreter','latex','FontSize',15,...
    'position',[-0.633319749858728,3.162301964601093,-1])
xlim([0 6])
ylim([10^-3 10^4])

%% Ploting all(WT7,WT10,ASL)

nbins = 1000;
data_HotWT7_1 = [duwrtlocalmaxmin_HotWT7{1}'/(VF_HotWT7.u_tau),...
    duNN_HotWT7{1}'/(VF_HotWT7.u_tau*res_WT7(1,1))]; % Example correlated non-Gaussian data

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

data_HotWT7_9 = [duwrtlocalmaxmin_HotWT7{9}'/(VF_HotWT7.u_tau),...
    duNN_HotWT7{9}'/(VF_HotWT7.u_tau*res_WT7(9,1))]; % Example correlated non-Gaussian data


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


data_HotWT10_2 = [duwrtlocalmaxmin_HotWT10{2}'/(VF_HotWT10.u_tau),...
    duNN_HotWT10{2}'/(VF_HotWT10.u_tau*res_WT10(2,1))]; % Example correlated non-Gaussian data

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

data_HotWT10_9 = [duwrtlocalmaxmin_HotWT10{9}'/(VF_HotWT10.u_tau),...
    duNN_HotWT10{9}'/(VF_HotWT10.u_tau*res_WT10(9,1))]; % Example correlated non-Gaussian data

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

data_SLPIVASL_7 = [duwrtlocalmaxmin_SLPIVASL{7}{25}'/(VF_SLPIVASL.u_tau),...
    duNN_SLPIVASL{7}{25}'/(VF_SLPIVASL.u_tau*res_ASL(7,25))]; % Example correlated non-Gaussian data

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


data_SLPIVASL_19 = [duwrtlocalmaxmin_SLPIVASL{19}{25}'/(VF_SLPIVASL.u_tau),...
    duNN_SLPIVASL{19}{25}'/(VF_SLPIVASL.u_tau*res_ASL(19,25))]; % Example correlated non-Gaussian data

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


%% 3D plot
nbins = 1000;
data_HotWT7_1 = [value_pairs_cell_HotWT7{1}(:,2)/VF_HotWT7.u_tau, value_pairs_cell_HotWT7{1}(:,1)/VF_HotWT7.u_tau]; % Example correlated non-Gaussian data
[~, X_edges, Y_edges] = histcounts2(data_HotWT7_1(:,1), data_HotWT7_1(:,2), nbins);

% Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = (Y_edges(1:end-1) + Y_edges(2:end)) / 2;
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);

% Plot 2D contour of joint PDF (flat on X-Y plane)
figure;
set(gcf,'Position',[743,355,806,519])
axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
scatter(data_HotWT7_1(:,1), data_HotWT7_1(:,2), 5, [0.0,1,0.0],'filled'); % Scatter original data
hold on;
% contourf(X_grid, Y_grid, pdf_hist_HotWT7_1', 100, 'LineStyle', 'none');
contour(X_grid, Y_grid, pdf_hist_HotWT7_1', 100, 'LineWidth', 3);
% Load and apply plasma colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
% cmap2 = [1,1,1;cmap2];

colormap(gca, cmap2);
caxis([0 0.5]);  % Match your original range

% Add colorbar (for j.p.d.f. values)
hcb = colorbar;
title(hcb, 'j.p.d.f.', 'Interpreter', 'Latex', 'FontSize', 20, 'position',...
    [-35.33734564345655,-3.758313275518872,0])
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
hcb.Ticks = [0,0.1,0.2,0.3,0.4,0.5];
hcb.TickLength = 0;
hcb.Location = "north";
hcb.AxisLocation = 'out';
hcb.FontSize = 14;
hcb.Position = [0.122828784119107,0.755298651252408,0.339950372208436,0.023121387283237];




% Set the range for conditional X
y_cond = 3.0;
y_lower = 2.95;
y_upper = 3.05;

% Extract conditional Y-values (X is in col 2, Y in col 1)
cond_indices = data_HotWT7_1(:,2) > y_lower & data_HotWT7_1(:,2) < y_upper;
cond_X = data_HotWT7_1(cond_indices, 1);



bin_edges = linspace(-6, 6, 35);  % You can change number of bins if needed
bin_width = bin_edges(2) - bin_edges(1);

% Compute histogram
counts = histcounts(cond_X, bin_edges);

% Normalize to get PDF (area under curve = 1)
x_pdf = counts / (sum(counts) * bin_width);

% Compute bin centers for plotting
x_vals = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% Compute conditional PDF of Y
% [x_pdf, x_vals] = ksdensity(cond_X,linspace(-6, 6, 1000));
% trapz(x_vals,x_pdf) or trapz(bin_width,x_pdf)
% Set the range for conditional X
x_cond = -1.0;
x_lower = -1.05;
x_upper = -.95;

% Extract conditional Y-values (X is in col 2, Y in col 1)
cond_indices = data_HotWT7_1(:,1) > x_lower & data_HotWT7_1(:,1) < x_upper;
cond_Y = data_HotWT7_1(cond_indices, 2);

bin_edges = linspace(-6, 6, 100);  % You can change number of bins if needed
bin_width = bin_edges(2) - bin_edges(1);

counts = histcounts(cond_Y, bin_edges);

% Normalize to get PDF (area under curve = 1)
y_pdf = counts / (sum(counts) * bin_width);

% Compute bin centers for plotting
y_vals = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
% Compute conditional PDF of Y
% [y_pdf, y_vals] = ksdensity(cond_Y,linspace(-6, 6, 1000));
% trapz(y_vals,y_pdf) or trapz(bin_width,y_pdf)

% Plot the conditional PDF in Z dimension
plot3(x_cond * ones(size(y_vals)),y_vals,  y_pdf, 'b-', 'LineWidth', 3);
plot3(x_vals,y_cond * ones(size(x_vals)),  x_pdf, 'r-', 'LineWidth', 3);


set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','ZGrid','on','XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on')
xlabel('$\frac{\mathrm{Local\:min}}{u_{\tau}}$', 'Interpreter', 'latex', 'FontSize', 22,'position',...
    [0.511648929807876,-8.238345813019066,-0.030586792187206],'Rotation',22);
ylabel('$\frac{\mathrm{Local\:max}}{u_{\tau}}$', 'Interpreter', 'latex', 'FontSize', 22,'position',...
    [-7.204871853750234,2,-0.091324770143437],'Rotation',314);
zlabel('p.d.f.', 'Interpreter', 'latex', 'FontSize', 16.5);
legend(sprintf('Pairwise local min-max $\\mathrm{u}^{\\prime}$ H-W(m1) z$/\\delta$ = %.2f',VF_HotWT7.z(1)/VF_HotWT7.delta),...
    'Histogram density estimation','$\mathrm{P}(\mathrm{min}\mid\mathrm{max1})$','$\mathrm{P}(\mathrm{max}\mid\mathrm{min1})$',...
    'Interpreter','latex','FontSize',11,'Position',...
    [0.007259247094543,0.847229537822918,0.443333881874271,0.146242769788915],...
    'Numcolumns',1,'Orientation','vertical','color','none');
annotation(gcf,'textbox',...
    [0.00596277915632756 0.643545279383429 0.0486277915632756 0.0616570327552985],...
    'String',{'(a)'},...
    'Interpreter','latex',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');
xlim([-6 6])
ylim([-6 6])
zlim([0 0.7])
axis square;
view(-31.245251406131388,36.259361702723517);



nbins = 1000;

data_HotWT7_1 = [duwrtlocalmaxmin_HotWT7{1}'/(VF_HotWT7.u_tau),...
    duNN_HotWT7{1}'/(VF_HotWT7.u_tau*res_WT7(1,1))]; % Example correlated non-Gaussian data

Y_min = min(data_HotWT7_1(:,2));
if Y_min <= 0
    Y_min = min(data_HotWT7_1(data_HotWT7_1(:,2) > 0, 2)); % Find smallest positive value
end

% Use `logspace` for Y bins (since we will use log scaling)
X_edges = linspace(min(data_HotWT7_1(:,1)), max(data_HotWT7_1(:,1)), nbins+1);
Y_edges = logspace(log10(Y_min), log10(max(data_HotWT7_1(:,2))), nbins+1);

% Compute 2D histogram
% [counts, X_edges, Y_edges] = histcounts2(data_HotWT7_1(:,1), data_HotWT7_1(:,2), X_edges, Y_edges);
[pdf_hist_HotWT7_1_normalized, X_edges, Y_edges] = histcounts2(data_HotWT7_1(:,1), data_HotWT7_1(:,2), X_edges, Y_edges,'Normalization','pdf');
%  Compute bin centers
X_centers = (X_edges(1:end-1) + X_edges(2:end)) / 2;
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[X_grid, Y_grid] = meshgrid(X_centers, Y_centers);


P_final = sum(pdf_hist_HotWT7_1_normalized .*(diff(X_edges)'* diff(Y_edges)),'all');
disp(['Normalized Total Probability (should be ≈ 1): ', num2str(P_final)]);


ax1 = axes('Position',[0.378347510861831,0.07097763082208,0.81791702174454,0.675282138848995]);
axes(ax1)
scatter(data_HotWT7_1(:,1), data_HotWT7_1(:,2), 5, [0.0,1,0.0],'filled'); % Scatter original data
hold on;
contour(X_grid, Y_grid, log(pdf_hist_HotWT7_1_normalized'), 600, 'LineWidth', 3); % Contour of PDF

% Apply custom colormap
cmap2 = getPyPlot_cMap('plasma', [], [],...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
cmap2 = flip(cmap2);
colormap(gca, cmap2);
caxis([-10 1])
hcb = colorbar;
title(hcb, 'log(j.p.d.f.)', 'Interpreter', 'Latex', 'FontSize', 20, 'position',...
    [-57.837345643456544,-3.758313275518872,0])
hcb.Label.Interpreter = 'latex';
hcb.TickLabelInterpreter = 'latex';
% hcb.Ticks = [0,0.2,0.4,0.6,0.8,1.0];
hcb.Ticks = [-10,-5,1];
hcb.TickLength = 0;
hcb.Location = "north";
hcb.AxisLocation = 'out';
hcb.FontSize = 14;
hcb.Position = [0.678660049627792,0.755298651252408,0.291563275434244,0.024259455915291];


x_cond = 0.1;
x_lower = x_cond-0.1;
x_upper = x_cond+0.1;

% Extract conditional Y-values (X is in col 2, Y in col 1)
cond_indices = data_HotWT7_1(:,1) > x_lower & data_HotWT7_1(:,1) < x_upper;
cond_Y = data_HotWT7_1(cond_indices, 2);


Y_edges = logspace(log10(1e-4), log10(1e4), 100);
Y_centers = sqrt(Y_edges(1:end-1) .* Y_edges(2:end)); % Use geometric mean for log bins
[y_pdf] = histcounts(cond_Y, Y_edges,'Normalization','pdf');
% histogram(cond_Y, Y_edges,'Normalization','count')

% dy = diff(Y_edges); 
% Normalize to get PDF (area under curve = 1)
% y_pdf = counts / (sum(counts(:)) * mean(dy));
% y_pdf = counts ./ (sum(counts .* diff(Y_edges)));

% Compute bin centers for plotting
y_vals = Y_centers;
% sum(y_pdf .* diff(Y_edges))


[counts_prob] = histcounts(cond_Y, Y_edges, 'Normalization', 'probability');
sum(counts_prob)
ax2 = axes('Position',[0.348555463545803,0.111440058567745,0.8028442914648,0.565455549253618]);
axes(ax2)
bar3(Y_centers,counts_prob)
set(gca,'YTickLabel',[],'TickLabelInterpreter','latex','FontSize',13,'Box','on','YScale','log',...
    'XGrid','on','YGrid','on','ZGrid','on','XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on','YDir', 'normal',...
    'color','none')
zlabel('Probability','Interpreter','latex','FontSize',13,'position',[-10.992832538325649,271082.4474212482,0.050000042037835]);
view(-20.485199999999999,21.699999999999999);

axes(ax1)
plot3(x_cond * ones(size(y_vals)), y_vals, y_pdf, 'b-', 'LineWidth', 3);



set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on','ZGrid','on','XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on','YScale', 'log',...
    'YTick', [10^-3,10^0,10^3])
xlabel('$\frac{\min\{\left|\mathrm{du}^{\prime}_{\mathrm{NLmin}}\right|, \left|\mathrm{du}^{\prime}_{\mathrm{NLmax}}\right|\}}{u_{\tau}}$',...
    'Interpreter','latex','FontSize',16,'position',[3.44777654669215,0.000199816178237,-0.001209669806895],'Rotation',20)
ylabel('$\frac{\left|\mathrm{du}^{\prime}_{\mathrm{NN}}\right|}{\Delta\mathrm{x}_{\mathrm{data}}u_{\tau}}$','Interpreter','latex','FontSize',22,...
    'position',[-1.132696750789691,6.584629037676163,-0.001528455674804],'Rotation',314);
zlabel('p.d.f.', 'Interpreter', 'latex', 'FontSize', 16.5);
text(-0.223849900087045,0.89982096746428,0,'(b)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',16.5);
legend(sprintf('$\\mathrm{Pairwise}\\ \\frac{\\min\\left\\{\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NLmin}}\\right|,\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NLmax}}\\right|\\right\\}}{u_{\\tau}},\\frac{\\left|\\mathrm{du}^{\\prime}_{\\mathrm{NN}}\\right|}{\\Delta\\mathrm{x}_{\\mathrm{data}}u_{\\tau}}$ H-W(m1) $\\mathrm{z}/\\delta=%.2f$', VF_HotWT7.z(1)/VF_HotWT7.delta),...
    'Histogram density estimation','$\mathrm{P}(\frac{\left|\mathrm{du}^{\prime}_{\mathrm{NN}}\right|}{\Delta\mathrm{x}_{\mathrm{data}}u_{\tau}}\mid\frac{\min\{\left|\mathrm{du}^{\prime}_{\mathrm{NLmin}}\right|, \left|\mathrm{du}^{\prime}_{\mathrm{NLmax}}\right|\}}{u_{\tau}})$',...
    'Interpreter','latex','FontSize',11.5,'Position',...
    [0.456206698508976,0.844209570953662,0.540605970975763,0.149137684811835],...
    'Numcolumns',1,'Orientation','vertical','color','none');

xlim([0 6])
ylim([10^-3 10^4])
% zlim([0 0.003])
axis square;
view(-31.245251406131388,36.259361702723517);
axes(ax2)


[counts_cdf] = histcounts(cond_Y, Y_edges, 'Normalization', 'cdf');
ax3 = axes('Position',[0.510702625819958,0.217413083615914,0.769174012768887,0.49801816967751]);
axes(ax3)
bar3(Y_centers,counts_cdf)
set(gca,'YTick',[1e-3,1e0,1e3],'TickLabelInterpreter','latex','FontSize',13,'Box','on','YScale','log',...
    'XGrid','on','YGrid','on','ZGrid','on','XMinorGrid','on','YMinorGrid','on','ZMinorGrid','on','YDir', 'normal',...
    'color',[1,1,1])
zlabel('c.d.f.','Interpreter','latex','FontSize',15);
ylim([1e-3 1e4])
view(-20.485199999999999,21.699999999999999);
annotation('arrow', [0.3, 0.5], [0.6, 0.8]);
text(1, 1, '$\int$', 'Interpreter', 'latex', 'Rotation', 90, 'FontSize', 12)

%%
Given_min = linspace(-5, 5, 100);  % Your input value as a scalar
mix_coeff_Given_min = zeros(size(Given_min,2),5);
mean_mix_Given_min = zeros(size(Given_min,2),5);
std_mix_Given_min = zeros(size(Given_min,2),5);
for i=1:size(Given_min,2)
    [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_min_model.py",...
        "ReturnList",x_test = Given_min(i));
    mix_coeff_Given_min(i,:) = result{1};
    mean_mix_Given_min(i,:) = result{2};
    std_mix_Given_min(i,:) = result{3};
end


Given_max = linspace(-5, 6, 120);  % Your input value as a scalar
mix_coeff_Given_max = zeros(size(Given_max,2),5);
mean_mix_Given_max = zeros(size(Given_max,2),5);
std_mix_Given_max = zeros(size(Given_max,2),5);

for i=1:size(Given_max,2)
    [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_max_model.py",...
        "ReturnList",x_test = Given_max(i))
    mix_coeff_Given_max(i,:) = result{1};
    mean_mix_Given_max(i,:) = result{2};
    std_mix_Given_max(i,:) = result{3};
end


Given_dist_NLmin_NLmax = linspace(0, 6, 100);
mix_coeff_Given_dist_NLmin_NLmax = zeros(size(Given_dist_NLmin_NLmax,2),5);
mean_mix_Given_dist_NLmin_NLmax = zeros(size(Given_dist_NLmin_NLmax,2),5);
std_mix_Given_dist_NLmin_NLmax = zeros(size(Given_dist_NLmin_NLmax,2),5);

for i=1:size(Given_dist_NLmin_NLmax,2)
    [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
        "ReturnList",x_test = Given_dist_NLmin_NLmax(i))
    mix_coeff_Given_dist_NLmin_NLmax(i,:) = result{1};
    mean_mix_Given_dist_NLmin_NLmax(i,:) = result{2};
    std_mix_Given_dist_NLmin_NLmax(i,:) = result{3};
end

figure
set(gcf,'Position',[823,262,806,960])
axes('Position',[0.066997518610422,0.66875,0.258064516129032,0.216666666666667])
plot(Given_min,mix_coeff_Given_min(:,1),'r-', 'LineWidth', 2)
hold on
plot(Given_min,mix_coeff_Given_min(:,2),'g--', 'LineWidth', 2)
plot(Given_min,mix_coeff_Given_min(:,3),'m:', 'LineWidth', 2)
plot(Given_min,mix_coeff_Given_min(:,4),'k-.', 'LineWidth', 2)
plot(Given_min,mix_coeff_Given_min(:,5),'b-', 'LineWidth', 1,'Marker','|',...
    'MarkerSize',7)
axis square
text(0.038461538461538,0.899038461538462,0,'(a)',...
    'Units','Normalized','Interpreter','LaTex','FontSize',16.5);
set(gca,'TickLabelInterpreter','latex','FontSize',14,'Box','on',...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on')
legend('1','2','3','4','5',...
    'Interpreter','latex','FontSize',12,'Position',...
    [0.291958457485421,0.968503244431069,0.366005115118567,0.024166666368643],...
    'Numcolumns',5,'Orientation','vertical','color','none');
xlabel('Local min sub--window u$^{\prime}/u_{\tau}$','Interpreter','latex','FontSize',14)
ylabel('$\mathrm{a}_{i}$','Interpreter','latex','FontSize',15)
xlim([-5 5])
ylim([0 0.5])


axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_min,mean_mix_Given_min(:,1))
hold on
plot(Given_min,mean_mix_Given_min(:,2))
plot(Given_min,mean_mix_Given_min(:,3))
plot(Given_min,mean_mix_Given_min(:,4))
plot(Given_min,mean_mix_Given_min(:,5))
axis square

axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_min,std_mix_Given_min(:,1))
hold on
plot(Given_min,std_mix_Given_min(:,2))
plot(Given_min,std_mix_Given_min(:,3))
plot(Given_min,std_mix_Given_min(:,4))
plot(Given_min,std_mix_Given_min(:,5))

axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_max,mix_coeff_Given_max(:,1))
hold on
plot(Given_max,mix_coeff_Given_max(:,2))
plot(Given_max,mix_coeff_Given_max(:,3))
plot(Given_max,mix_coeff_Given_max(:,4))
plot(Given_max,mix_coeff_Given_max(:,5))

axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_max,mean_mix_Given_max(:,1))
hold on
plot(Given_max,mean_mix_Given_max(:,2))
plot(Given_max,mean_mix_Given_max(:,3))
plot(Given_max,mean_mix_Given_max(:,4))
plot(Given_max,mean_mix_Given_max(:,5))

axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_max,std_mix_Given_max(:,1))
hold on
plot(Given_max,std_mix_Given_max(:,2))
plot(Given_max,std_mix_Given_max(:,3))
plot(Given_max,std_mix_Given_max(:,4))
plot(Given_max,std_mix_Given_max(:,5))

figure
% axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_dist_NLmin_NLmax,mix_coeff_Given_dist_NLmin_NLmax(:,1))
hold on
plot(Given_dist_NLmin_NLmax,mix_coeff_Given_dist_NLmin_NLmax(:,2))
plot(Given_dist_NLmin_NLmax,mix_coeff_Given_dist_NLmin_NLmax(:,3))
plot(Given_dist_NLmin_NLmax,mix_coeff_Given_dist_NLmin_NLmax(:,4))
plot(Given_dist_NLmin_NLmax,mix_coeff_Given_dist_NLmin_NLmax(:,5))

figure
% axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_dist_NLmin_NLmax,mean_mix_Given_dist_NLmin_NLmax(:,1))
hold on
plot(Given_dist_NLmin_NLmax,mean_mix_Given_dist_NLmin_NLmax(:,2))
plot(Given_dist_NLmin_NLmax,mean_mix_Given_dist_NLmin_NLmax(:,3))
plot(Given_dist_NLmin_NLmax,mean_mix_Given_dist_NLmin_NLmax(:,4))
plot(Given_dist_NLmin_NLmax,mean_mix_Given_dist_NLmin_NLmax(:,5))

figure
% axes('Position',[-0.125374573505415,0.07097763082208,0.81791702174454,0.675282138848995])
plot(Given_dist_NLmin_NLmax,std_mix_Given_dist_NLmin_NLmax(:,1))
hold on
plot(Given_dist_NLmin_NLmax,std_mix_Given_dist_NLmin_NLmax(:,2))
plot(Given_dist_NLmin_NLmax,std_mix_Given_dist_NLmin_NLmax(:,3))
plot(Given_dist_NLmin_NLmax,std_mix_Given_dist_NLmin_NLmax(:,4))
plot(Given_dist_NLmin_NLmax,std_mix_Given_dist_NLmin_NLmax(:,5))

