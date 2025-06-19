classdef VFfeaturedVorX_mac < dynamicprops
    properties
        name, type, u, uprime, w, x, z, delta, u_tau
        Lambda_ci, LSM, Delx, HRVFu, HRVFw, HRVFx
        HRVFz, factor, HRVFuprime, Rij, omega, Lambda_cirms, z0%dudx, dudz
        %dwdx, dwdz
    end
    
    methods
        function Obj = VFfeaturedVorX_mac(name, type, u, uprime, w, x, z, delta, u_tau,...
                LSM, Delx, Lambda_cirms, z0)
            
            Obj.name = name;
            Obj.type = type;
            Obj.u = u;
            Obj.uprime = uprime;
            Obj.w = w;
            Obj.x = x;
            Obj.z = z;
            Obj.delta = delta;
            Obj.u_tau = u_tau;
            Obj.LSM = LSM;
            Obj.Delx = Delx;
            Obj.Lambda_cirms = Lambda_cirms;
            Obj.z0 = z0;
        end
        
        function Lambdaci(Obj)
            if Obj.type == "Exp"
                Obj.Lambda_ci = zeros(size(Obj.u));
                %                 Obj.dudx = zeros(size(Obj.u));
                %                 Obj.dudz = zeros(size(Obj.u));
                %                 Obj.dwdx = zeros(size(Obj.u));
                %                 Obj.dwdz = zeros(size(Obj.u));
                DelxExp = Obj.x(2)-Obj.x(1);
                DelzExp = Obj.z(2)-Obj.z(1);
            else
                %                 L = ceil(Obj.LSM*Obj.delta/Obj.Delx);%Number of intervals
                %                 Numrator = floor((size(Obj.u,2)-1)/L);
                %                 Obj.Lambda_ci = zeros(size(Obj.u,1),L+1,Numrator);
                %                 Obj.dudx = zeros(size(Obj.u,1),L+1,Numrator);
                %                 Obj.dudz = zeros(size(Obj.u,1),L+1,Numrator);
                %                 Obj.dwdx = zeros(size(Obj.u,1),L+1,Numrator);
                %                 Obj.dwdz = zeros(size(Obj.u,1),L+1,Numrator);
                %                 DelxGen = Obj.Delx;
                %                 DelzGen = Obj.z(2)-Obj.z(1);
                %                 Obj.dudx = zeros(size(Obj.HRVFu));
                %                 Obj.dudz = zeros(size(Obj.HRVFu));
                %                 Obj.dwdx = zeros(size(Obj.HRVFu));
                %                 Obj.dwdz = zeros(size(Obj.HRVFu));
                if isempty(Obj.Lambda_ci) 
                    Obj.Lambda_ci = zeros(size(Obj.HRVFu));
                    Obj.omega = zeros(size(Obj.HRVFu));
                end
                DelxGen = Obj.HRVFx(2)-Obj.HRVFx(1);
                DelzGen = Obj.z(2)-Obj.z(1);
            end
            if Obj.type == "Exp"
                for S= 1:size(Obj.Lambda_ci, 3)
                    
                    %                     for r = 2:size(Obj.Lambda_ci, 1)-1
                    %
                    %                         for c = 2:size(Obj.Lambda_ci, 2)-1
                    %
                    % %                             Obj.dudx(r,c,S) = (Obj.u(r, c+1, S)-Obj.u(r, c-1, S))/(2*DelxExp);
                    % %                             Obj.dudz(r,c,S) = (Obj.u(r+1, c, S)-Obj.u(r-1, c, S))/(2*DelzExp);
                    % %                             Obj.dwdx(r,c,S) = (Obj.w(r, c+1, S)-Obj.w(r, c-1, S))/(2*DelxExp);
                    % %                             Obj.dwdz(r,c,S) = (Obj.w(r+1, c, S)-Obj.w(r-1, c, S))/(2*DelzExp);
                    % %                             Mat=[Obj.dudx(r,c,S) Obj.dudz(r,c,S);...
                    % %                                 Obj.dwdx(r,c,S) Obj.dwdz(r,c,S)];
                    % %                             Obj.omega(r,c,S) = Obj.dwdx(r,c,S) - Obj.dudz(r,c,S);
                    %                             dudx = (Obj.u(r, c+1, S)-Obj.u(r, c-1, S))/(2*DelxExp);
                    %                             dudz = (Obj.u(r+1, c, S)-Obj.u(r-1, c, S))/(2*DelzExp);
                    %                             dwdx = (Obj.w(r, c+1, S)-Obj.w(r, c-1, S))/(2*DelxExp);
                    %                             dwdz = (Obj.w(r+1, c, S)-Obj.w(r-1, c, S))/(2*DelzExp);
                    %                             Mat=[dudx dudz;...
                    %                                 dwdx dwdz];
                    %                             Obj.omega(r,c,S) = dwdx - dudz;
                    %                             Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(Mat))))...
                    %                                 *sign(Obj.omega(r,c,S));
                    %                         end
                    %                     end
                    
                    [dudx, dudz] = gradient(Obj.u(:,:, S),DelxExp,DelzExp);
                    [dwdx, dwdz] = gradient(Obj.w(:,:, S),DelxExp,DelzExp);
                    Obj.omega(:,:,S) = dwdx - dudz;
                    Mat = zeros(size(dudx,1),size(dudx,2),2,2);
                    Mat(:,:,1,1) = dudx;
                    Mat(:,:,1,2) = dudz;
                    Mat(:,:,2,1) = dwdx;
                    Mat(:,:,2,2) = dwdz;
                    for r = 1:size(Obj.Lambda_ci,1)
                        
                        for c = 1:size(Obj.Lambda_ci,2)
                            temp = reshape(Mat(r,c,:,:),[2,2]);
                            Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                                *sign(Obj.omega(r,c,S));
                            
                        end
                    end
                    
                end
                Obj.Lambda_ci(find(abs(Obj.Lambda_ci)<Obj.Lambda_cirms))=0;
            else
                for S= 1:size(Obj.Lambda_ci,3)
                    
                    %                     ucroped = Obj.u(:, (S-1)*L+1:(S)*L+1);
                    %                     wcroped = Obj.w(:, (S-1)*L+1:(S)*L+1);
                    
                    %                     for r = 2:size(Obj.Lambda_ci,1)-1
                    %
                    %                         for c = 2:size(Obj.Lambda_ci,2)-1
                    %
                    % %                             Obj.dudx(r,c,S) = (ucroped(r, c+1)- ucroped(r, c-1))/(2*DelxGen);
                    % %                             Obj.dudz(r,c,S) = (ucroped(r+1, c)- ucroped(r-1, c))/(2*DelzGen);
                    % %                             Obj.dwdx(r,c,S) = (wcroped(r, c+1)- wcroped(r, c-1))/(2*DelxGen);
                    % %                             Obj.dwdz(r,c,S) = (wcroped(r+1, c)- wcroped(r-1, c))/(2*DelzGen);
                    %
                    % %                             Obj.dudx(r,c,S) = (Obj.HRVFu(r, c+1, S)- Obj.HRVFu(r, c-1, S))/(2*DelxGen);
                    % %                             Obj.dudz(r,c,S) = (Obj.HRVFu(r+1, c, S)- Obj.HRVFu(r-1, c, S))/(2*DelzGen);
                    % %                             Obj.dwdx(r,c,S) = (Obj.HRVFw(r, c+1, S)- Obj.HRVFw(r, c-1, S))/(2*DelxGen);
                    % %                             Obj.dwdz(r,c,S) = (Obj.HRVFw(r+1, c, S)- Obj.HRVFw(r-1, c, S))/(2*DelzGen);
                    % %                             Mat=[Obj.dudx(r,c,S) Obj.dudz(r,c,S);...
                    % %                                 Obj.dwdx(r,c,S) Obj.dwdz(r,c,S)];
                    % %                             Obj.omega(r,c,S) = Obj.dwdx(r,c,S) - Obj.dudz(r,c,S);
                    %                             dudx = (Obj.HRVFu(r, c+1, S)- Obj.HRVFu(r, c-1, S))/(2*DelxGen);
                    %                             dudz = (Obj.HRVFu(r+1, c, S)- Obj.HRVFu(r-1, c, S))/(2*DelzGen);
                    %                             dwdx = (Obj.HRVFw(r, c+1, S)- Obj.HRVFw(r, c-1, S))/(2*DelxGen);
                    %                             dwdz = (Obj.HRVFw(r+1, c, S)- Obj.HRVFw(r-1, c, S))/(2*DelzGen);
                    %                             Mat=[dudx dudz;...
                    %                                 dwdx dwdz];
                    %                             Obj.omega(r,c,S) = dwdx - dudz;
                    %                             Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(Mat))))...
                    %                                 *sign(Obj.omega(r,c,S));
                    %
                    %                         end
                    %                     end
                    
                    [dudx, dudz] = gradient(Obj.HRVFu(:,:, S),DelxGen,DelzGen);
                    [dwdx, dwdz] = gradient(Obj.HRVFw(:,:, S),DelxGen,DelzGen);
                    Obj.omega(:,:,S) = dwdx - dudz;
                    Mat = zeros(size(dudx,1),size(dudx,2),2,2);
                    Mat(:,:,1,1) = dudx;
                    Mat(:,:,1,2) = dudz;
                    Mat(:,:,2,1) = dwdx;
                    Mat(:,:,2,2) = dwdz;
                    for r = 1:size(Obj.Lambda_ci,1)
                        
                        for c = 1:size(Obj.Lambda_ci,2)
                            temp = reshape(Mat(r,c,:,:),[2,2]);
                            Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                                *sign(Obj.omega(r,c,S));
                            
                        end
                    end
                    
                    
                end
                Obj.Lambda_ci(find(abs(Obj.Lambda_ci)<Obj.Lambda_cirms))=0;
            end
            
        end
        function Increasing_Res(Obj, factor)
            if Obj.type == "Exp"
                
            else
                L = ceil(Obj.LSM*Obj.delta/Obj.Delx);%Number of intervals
                Numrator = floor((size(Obj.u,2)-1)/L);
                DelxGen = Obj.Delx;
                DelzGen = Obj.z(2)-Obj.z(1);
                if isnan(factor)
                    Obj.factor = floor(DelxGen/DelzGen);
                else
                    Obj.factor = factor;
                end
                Obj.HRVFu=zeros(size(Obj.u,1), Obj.factor*(L) +1, Numrator);
                Obj.HRVFw=zeros(size(Obj.u,1), Obj.factor*(L) +1, Numrator);
                Obj.HRVFx= 0:DelxGen/Obj.factor: (L)*DelxGen;
                Obj.HRVFz=Obj.z;
                xorg = 0:DelxGen:(L)*DelxGen;
                Obj.x = xorg;
                xi = Obj.HRVFx;
                kappa = 0.39;
                for S= 1:Numrator
                    
                    ucroped = Obj.u(:, (S-1)*L+1:(S)*L+1);
                    wcroped = Obj.w(:, (S-1)*L+1:(S)*L+1);
                    
                    for r = 1:size(Obj.HRVFz,1)
                        % Interpolate the row
                        Obj.HRVFu(r, :, S) = interp1(xorg, ucroped(r, :), xi, 'linear');
                        Obj.HRVFw(r, :, S) = interp1(xorg, wcroped(r, :), xi, 'linear');
                    end
                end
%                 Obj.HRVFuprime = Obj.HRVFu-mean(Obj.HRVFu,3);
                Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);
            end
            
        end
        function Gauss_kernel(Obj,kern_xsize, kern_ysize)
            mu = [0 0];
            Sigma = [Obj.Delx/5000 Obj.Delx/5000*sqrt(tand(12.9));...
                Obj.Delx/5000*sqrt(tand(12.9)) tand(13)*Obj.Delx/5000];
            %             Sigma = [Obj.Delx/100 sqrt(Obj.Delx/100*Obj.Delx/100*tand(0));...
            %                 sqrt(Obj.Delx/100*Obj.Delx/100*tand(0)) Obj.Delx/100];
            [X,Y] = meshgrid(-kern_xsize/2*Obj.Delx:Obj.Delx/Obj.factor:kern_xsize/2*Obj.Delx,...
                -kern_ysize/2*Obj.Delx:Obj.Delx/Obj.factor:kern_ysize/2*Obj.Delx);
            kernel = mvnpdf([X(:) Y(:)],mu,Sigma);
            kernel = reshape(kernel,size(X));
            
            figure
            set(gcf,'Position',[806,788,560,271])
            contourf(X/Obj.Delx, Y/Obj.Delx, kernel, 100,'LineStyle','none')
            colorbar
            hcb2=colorbar;
            hcb2.TickLabelInterpreter = 'latex';
            set(gca,'TickLabelInterpreter','latex','FontSize',13)
            xlabel('x/$\lambda_{T}$','Interpreter','Latex','FontSize',14);
            ylabel('z/$\lambda_{T}$','Interpreter','Latex','FontSize',14);
            axis equal
            %check the 2D-integral of the p should be 1 for very small std
            %                 integralOverX = trapz(X(1, :), p, 2); % Integration over columns
            %
            %                 % Integrate the result over y-direction
            %                 totalIntegral = trapz(Y(:, 1), integralOverX);
            %
            %                 % Display the result
            %                 disp(['Integral of p over the grid: ', num2str(totalIntegral)]);
            normalization_matrix = conv2(ones(size(Obj.HRVFu,1),size(Obj.HRVFu,2)),...
                kernel, 'same');
            kappa = 0.39;
            for S =1:size(Obj.HRVFu,3)
                
                Obj.HRVFu(:,:,S) = conv2(Obj.HRVFu(:,:,S), kernel, 'same')./normalization_matrix;
                Obj.HRVFw(:,:,S) = conv2(Obj.HRVFw(:,:,S), kernel, 'same')./normalization_matrix;
%                 Obj.HRVFuprime(:,:,S) = conv2(Obj.HRVFuprime(:,:,S), kernel, 'same')./normalization_matrix;
            end
%             Obj.HRVFuprime = Obj.HRVFu-mean(Obj.HRVFu,3);
            Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);
        end
        function Twop_corr(Obj,zref)
            if Obj.type == "Exp"
                
                Obj.Rij = zeros(size(Obj.u,1),2*size(Obj.u,2)-1,size(Obj.u,3));
                Ind = find(Obj.z/Obj.delta > zref, 1);
                stduref = mean(rms(Obj.uprime(Ind,:,:),3),2);
                stdutargetvec = mean(rms(Obj.uprime(:,:,:),3),2);
                for j=1:size(Obj.u,3)
                    
                    uprimeref = Obj.uprime(Ind,:,j);
                    
                    
                    for i=1:size(Obj.u,1)
                        
                        uprimetarget = Obj.uprime(i,:,j);
                        stdutarget = stdutargetvec(i,1);
                        [r,~]=xcorr(uprimeref,uprimetarget,'unbiased');
                        Obj.Rij(i,:,j) = r./(stduref*stdutarget);
                        
                    end
                end
                
            else
                
                Obj.Rij = zeros(size(Obj.HRVFu,1),2*size(Obj.HRVFu,2)-1,size(Obj.HRVFu,3));
                Ind = find(Obj.z/Obj.delta > zref, 1);
                stduref = mean(rms(Obj.HRVFuprime(Ind,:,:),3),2);
                stdutargetvec = mean(rms(Obj.HRVFuprime(:,:,:),3),2);
                %                 stduref = rms(reshape(Obj.HRVFuprime(Ind,:,:),...
                %                     [1,size(Obj.HRVFuprime,2)*size(Obj.HRVFuprime,3)]),2);
                %                 stdutargetvec = rms(reshape(Obj.HRVFuprime(:,:,:),...
                %                     [size(Obj.HRVFuprime,1),size(Obj.HRVFuprime,2)*size(Obj.HRVFuprime,3)]),2);
                for S =1:size(Obj.HRVFu,3)
                    
                    uprimeref = Obj.HRVFuprime(Ind,:,S);
                    
                    for i=1:size(Obj.z,1)
                        
                        uprimetarget = Obj.HRVFuprime(i,:,S);
                        stdutarget = stdutargetvec(i,1);
                        [r,~]=xcorr(uprimeref,uprimetarget,'unbiased');
                        Obj.Rij(i,:,S) = r./(stduref*stdutarget);
                        
                    end
                end
            end
        end
        function VortXmodeling(Obj)
            if Obj.name == "GenWT7"
                WT7VorX = load(['/Users/roozbehehsani/Google Drive/My Drive/Research/DATABASE/' ...
                    'Vortex statistics/data/mesh_07ms_vortex_properties.mat']);
                WT7VorX_scales = load(['/Users/roozbehehsani/Google Drive/My Drive/Research/DATABASE/' ...
                    'Vortex statistics/data/mesh_07ms_scales.mat']);
                R1 = corrcoef(WT7VorX.d_pro,WT7VorX.dU_pro);
                %                 R2 = corrcoef(WT7VorX.d_ret,WT7VorX.dU_ret);
                
                differencespro = abs(WT7VorX.z_pro - WT7VorX_scales.z.');
                [~, resultpro] = min(differencespro, [], 2);
                
                %                 differencesret = abs(WT7VorX.z_ret - WT7VorX_scales.z.');
                %                 [~, resultret] = min(differencesret, [], 2);
                %
                %                 figure
                %                 subplot(1,2,1)
                %                 plot(WT7VorX.d_pro./WT7VorX_scales.lambda(resultpro),...
                %                     -2*WT7VorX.dU_pro/Obj.u_tau,'b.')
                %                 set(gca,'TickLabelInterpreter','latex','FontSize',13)
                %                 xlabel('$d_{w}/\lambda_{T}$','Interpreter','Latex','FontSize',14);
                %                 ylabel('$\Delta u_{w}/u_{\tau}$','Interpreter','Latex','FontSize',14);
                %                 legend('Prograde','Interpreter','latex','FontSize',9)
                %                 xlim([0 4])
                %                 ylim([0 11])
                %                 subplot(1,2,2)
                %                 plot(WT7VorX.d_ret./WT7VorX_scales.lambda(resultret),...
                %                     2*WT7VorX.dU_ret/Obj.u_tau,'r.')
                %                 set(gca,'TickLabelInterpreter','latex','FontSize',13)
                %                 xlabel('$d_{w}/\lambda_{T}$','Interpreter','Latex','FontSize',14);
                %                 ylabel('$\Delta u_{w}/u_{\tau}$','Interpreter','Latex','FontSize',14);
                %                 legend('Retrograde','Interpreter','latex','FontSize',9)
                %                 xlim([0 4])
                %                 ylim([0 11])
                %
                
                
                %                 [pdfd_proOlam,d_proOlamp] = ksdensity(WT7VorX.d_pro./WT7VorX_scales.lambda(resultpro),...
                %                     'NumPoints',50);
                %                 [pdfd_retOlam,d_retOlamp] = ksdensity(WT7VorX.d_ret./WT7VorX_scales.lambda(resultret),...
                %                     'NumPoints',50);
                %
                %                 [cdfd_proOlam,d_proOlamc] = ksdensity(WT7VorX.d_pro./WT7VorX_scales.lambda(resultpro),...
                %                     'Function', 'cdf','NumPoints',50);
                %                 [cdfd_retOlam,d_retOlamc] = ksdensity(WT7VorX.d_ret./WT7VorX_scales.lambda(resultret),...
                %                     'Function', 'cdf','NumPoints',50);
                %
                %
                %
                %                 [pdfdU_proOutau,dU_proOutaup] = ksdensity(-2*WT7VorX.dU_pro/Obj.u_tau,...
                %                     'NumPoints',50);
                %                 [pdfdU_retOutau,dU_retOutaup] = ksdensity(2*WT7VorX.dU_ret/Obj.u_tau,...
                %                     'NumPoints',50);
                %
                %                 [cdfdU_proOutau,dU_proOutauc] = ksdensity(-2*WT7VorX.dU_pro/Obj.u_tau,...
                %                     'Function', 'cdf','NumPoints',50);
                %                 [cdfdU_retOutau,dU_retOutauc] = ksdensity(2*WT7VorX.dU_ret/Obj.u_tau,...
                %                     'Function', 'cdf','NumPoints',50);
                
                % Find Near-wall vortices
                m = 2;
                Near_wall_d_pro = WT7VorX.d_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*30*Obj.z0);
                Near_wall_dU_pro = WT7VorX.dU_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*30*Obj.z0);
                Near_wall_d_ret = WT7VorX.d_ret (WT7VorX.z_ret>=0 & WT7VorX.z_ret <= m*30*Obj.z0);
                Near_wall_dU_ret = WT7VorX.dU_ret (WT7VorX.z_ret>=0 & WT7VorX.z_ret <= m*30*Obj.z0);
                
                Near_wall_diffpro = abs(WT7VorX.z_pro(WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*30*Obj.z0)...
                    - WT7VorX_scales.z.');
                [~, result_Near_wall_pro] = min(Near_wall_diffpro, [], 2);
                
                Near_wall_diffret = abs(WT7VorX.z_ret(WT7VorX.z_ret>=0 & WT7VorX.z_ret <= m*30*Obj.z0)...
                    - WT7VorX_scales.z.');
                [~, result_Near_wall_ret] = min(Near_wall_diffret, [], 2);
                
                [pdfd_proOlam_Near_wall,d_proOlamp_Near_wall] = ksdensity(Near_wall_d_pro./...
                    WT7VorX_scales.lambda(result_Near_wall_pro),...
                    'NumPoints',50);
                [pdfd_retOlam_Near_wall,d_retOlamp_Near_wall] = ksdensity(Near_wall_d_ret./...
                    WT7VorX_scales.lambda(result_Near_wall_ret),...
                    'NumPoints',50);
                
                [cdfd_proOlam_Near_wall,d_proOlamc_Near_wall] = ksdensity(Near_wall_d_pro./...
                    WT7VorX_scales.lambda(result_Near_wall_pro),...
                    'Function', 'cdf','NumPoints',50);
                [cdfd_retOlam_Near_wall,d_retOlamc_Near_wall] = ksdensity(Near_wall_d_ret./...
                    WT7VorX_scales.lambda(result_Near_wall_ret),...
                    'Function', 'cdf','NumPoints',50);
                
                
                
                [pdfdU_proOutau_Near_wall,dU_proOutaup_Near_wall] = ksdensity(-2*Near_wall_dU_pro/Obj.u_tau,...
                    'NumPoints',50);
                [pdfdU_retOutau_Near_wall,dU_retOutaup_Near_wall] = ksdensity(2*Near_wall_dU_ret/Obj.u_tau,...
                    'NumPoints',50);
                
                [cdfdU_proOutau_Near_wall,dU_proOutauc_Near_wall] = ksdensity(-2*Near_wall_dU_pro/Obj.u_tau,...
                    'Function', 'cdf','NumPoints',50);
                [cdfdU_retOutau_Near_wall,dU_retOutauc_Near_wall] = ksdensity(2*Near_wall_dU_ret/Obj.u_tau,...
                    'Function', 'cdf','NumPoints',50);
                
                
            elseif Obj.name == "GenWT10"
                WT10VorX = load(['/Users/roozbehehsani/Google Drive/My Drive/Research/DATABASE/...' ...
                    'Vortex statistics/data/mesh_10ms_vortex_properties.mat']);
                R1 = corrcoef(WT10VorX.d_pro,WT10VorX.dU_pro);
%                 R2 = corrcoef(WT10VorX.d_ret,WT10VorX.dU_ret);
            elseif Obj.name == "GenASL"
                ASLVorX = load(['/Users/roozbehehsani/Google Drive/My Drive/Research/DATABASE/...' ...
                    'Vortex statistics/data/asl_vortex_properties.mat']);
                R1 = corrcoef(ASLVorX.d_pro,ASLVorX.dU_pro);
%                 R2 = corrcoef(ASLVorX.d_ret,ASLVorX.dU_ret);
            end
            kappa = 0.39;
            rng(123); % Seed for random number generator
            r = abs(R1(1,2)); % Target (Spearman) correlation %-0.4
            n = 10000000; % Number of samples
            
            % Functions
            %             gen.gauss.cop = @(r, n) genGaussCop(r, n);
            
            % Data generation and visualization
            randomvals = genGaussCop(r, n);
            rdelUw = randomvals(:,1);
            rdw = randomvals(:,2);
            ii = 1;
            C1 = trapz((1:0.1:9),exp(-(1:0.1:9)./0.85));
            %             trapz((1:0.1:9),1/C1*exp(-(1:0.1:9)./0.85))
            C2 = 1/(C1*0.85*(exp(-1/0.85)-exp(-9/0.85)));
            delUwOu_tau = -0.85*log(-rdelUw./(C2*C1*0.85)+exp(-1/0.85));
            
            
            
            %             figure
            %             subplot(1,2,1)
            %             plot(dU_proOutaup,pdfdU_proOutau,...
            %                 'LineStyle','none','color','b','Marker','^');
            %             hold on
            %             plot(dU_retOutaup,pdfdU_retOutau,...
            %                 'LineStyle','none','color','r','Marker','*');
            %             semilogy((1:0.1:9),1/C1*exp(-(1:0.1:9)./0.85),'k','LineWidth',1.5)
            %             set(gca,'YScale','log','TickLabelInterpreter','latex','FontSize',13)
            %             xlabel('$\Delta u_{w}/u_{\tau}$','Interpreter','Latex','FontSize',14);
            %             ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
            %             legend('Prograge','Retrograde','Gen','Interpreter','latex','FontSize',9)
            %             ylim([1e-5 1e1])
            %             subplot(1,2,2)
            %             plot(dU_proOutauc,cdfdU_proOutau,...
            %                 'LineStyle','none','color','b','Marker','^');
            %             hold on
            %             plot(dU_retOutauc,cdfdU_retOutau,...
            %                 'LineStyle','none','color','r','Marker','*');
            %             plot((1:0.001:9),-C2*C1*0.85*(-exp(-1/0.85)+exp(-(1:0.001:9)./0.85)),...
            %                 'k','LineWidth',1.5)
            % %             plot(delUwOu_tau,rdelUw,'Y.')
            %             set(gca,'TickLabelInterpreter','latex','FontSize',13)
            %             xlabel('$\Delta u_{w}/u_{\tau}$','Interpreter','Latex','FontSize',14);
            %             ylabel('c.d.f.','Interpreter','Latex','FontSize',14);
            %             legend('Prograge','Retrograde','Gen','Interpreter','latex','FontSize',9)
            
            
            data = double(WT7VorX.d_pro./WT7VorX_scales.lambda(resultpro));
            data = data(data>0.1);
            
            
            
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
            lb = [mean(log(data)), 0.01, 5, prctile(data, 75)];
            ub = [0, Inf, 6, prctile(data, 90)];
            
            % Optimize using fmincon
            options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
            [opt_params, ~] = fmincon(objective_log_pareto, initial_params, [], [], [], [], lb, ub, [], options);
            
            % Extract optimized parameters
            mu_opt = opt_params(1);
            sigma_opt = opt_params(2);
            alpha_opt = opt_params(3);
            x_t_opt = opt_params(4);
            
            % Display the results
            %             disp(['Optimized mu: ', num2str(mu_opt)]);
            %             disp(['Optimized sigma: ', num2str(sigma_opt)]);
            %             disp(['Optimized alpha: ', num2str(alpha_opt)]);
            %             disp(['Optimized x_t: ', num2str(x_t_opt)]);
            
            x_values = linspace(min(data), max(data), 100);
            fitted_pdf = log_pareto_pdf(x_values, mu_opt, sigma_opt,...
                alpha_opt,x_t_opt);
            fitted_cdf = log_pareto_cdf(x_values, mu_opt, sigma_opt,...
                alpha_opt,x_t_opt);
            
            
            
            F_xt = normcdf((log(x_t_opt) - mu_opt) / sigma_opt); % Log-normal CDF at x_t
            dwOlambda_T = zeros(size(rdw));
            
            % For u <= F_xt (Log-normal region)
            idx_log_normal = rdw <= F_xt;
            dwOlambda_T(idx_log_normal) = exp(mu_opt + sigma_opt * norminv(rdw(idx_log_normal)));
            
            % For u > F_xt (Power-law region)
            idx_power_law = rdw > F_xt;
            dwOlambda_T(idx_power_law) = x_t_opt * (1 - (rdw(idx_power_law) - F_xt) / (1 - F_xt)).^(-1 / alpha_opt);
            
            %             figure
            %             subplot(1,2,1)
            %             plot(d_proOlamp,pdfd_proOlam,...
            %                 'LineStyle','none','color','b','Marker','^');
            %             hold on
            %             plot(d_retOlamp,pdfd_retOlam,...
            %                 'LineStyle','none','color','r','Marker','*');
            %             loglog(x_values,fitted_pdf,'k','LineWidth',1.5)
            %             set(gca,'Xscale','log','YScale','log','TickLabelInterpreter','latex','FontSize',13)
            %             xlabel('$d_{w}/\lambda_{T}$','Interpreter','Latex','FontSize',14);
            %             ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
            %             legend('Prograge','Retrograde','Gen','Interpreter','latex','FontSize',9)
            %             ylim([1e-5 1e1])
            %             subplot(1,2,2)
            %             plot(d_proOlamc,cdfd_proOlam,...
            %                 'LineStyle','none','color','b','Marker','^');
            %             hold on
            %             plot(d_retOlamc,cdfd_retOlam,...
            %                 'LineStyle','none','color','r','Marker','*');
            %             plot(x_values,fitted_cdf,...
            %                 'k','LineWidth',1.5)
            % %             plot(dwOlambda_T,rdw,'y.')
            %             set(gca,'TickLabelInterpreter','latex','FontSize',13)
            %             xlabel('$d_{w}/\lambda_{T}$','Interpreter','Latex','FontSize',14);
            %             ylabel('c.d.f.','Interpreter','Latex','FontSize',14);
            %             legend('Prograge','Retrograde','Gen','Interpreter','latex','FontSize',9)
            desired_corr = -0.4;
            Trans = [1, desired_corr; 0, sqrt(1 - desired_corr^2)];
            
            % VorX near the wall
            figure
            subplot(1,2,1)
            plot(dU_proOutaup_Near_wall,pdfdU_proOutau_Near_wall,...
                'LineStyle','none','color','b','Marker','^');
            hold on
            plot(dU_retOutaup_Near_wall,pdfdU_retOutau_Near_wall,...
                'LineStyle','none','color','r','Marker','*');
            semilogy((1:0.1:9),1/C1*exp(-(1:0.1:9)./0.85),'k','LineWidth',1.5)
            set(gca,'YScale','log','TickLabelInterpreter','latex','FontSize',13)
            xlabel('$\Delta u_{w}/u_{\tau}$','Interpreter','Latex','FontSize',14);
            ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
            legend('Near wall Prograge','Near wall Retrograde','Gen','Interpreter','latex','FontSize',9)
            ylim([1e-5 1e1])
            subplot(1,2,2)
            plot(dU_proOutauc_Near_wall,cdfdU_proOutau_Near_wall,...
                'LineStyle','none','color','b','Marker','^');
            hold on
            plot(dU_retOutauc_Near_wall,cdfdU_retOutau_Near_wall,...
                'LineStyle','none','color','r','Marker','*');
            plot((1:0.001:9),-C2*C1*0.85*(-exp(-1/0.85)+exp(-(1:0.001:9)./0.85)),...
                'k','LineWidth',1.5)
            set(gca,'TickLabelInterpreter','latex','FontSize',13)
            xlabel('$\Delta u_{w}/u_{\tau}$','Interpreter','Latex','FontSize',14);
            ylabel('c.d.f.','Interpreter','Latex','FontSize',14);
            legend('Near wall Prograge','Near wall Retrograde','Gen','Interpreter','latex','FontSize',9)
            
            figure
            subplot(1,2,1)
            plot(d_proOlamp_Near_wall,pdfd_proOlam_Near_wall,...
                'LineStyle','none','color','b','Marker','^');
            hold on
            plot(d_retOlamp_Near_wall,pdfd_retOlam_Near_wall,...
                'LineStyle','none','color','r','Marker','*');
            loglog(x_values,fitted_pdf,'k','LineWidth',1.5)
            set(gca,'Xscale','log','YScale','log','TickLabelInterpreter','latex','FontSize',13)
            xlabel('$d_{w}/\lambda_{T}$','Interpreter','Latex','FontSize',14);
            ylabel('p.d.f.','Interpreter','Latex','FontSize',14);
            legend('Near wall Prograge','Near wall Retrograde','Gen','Interpreter','latex','FontSize',9)
            ylim([1e-5 1e1])
            subplot(1,2,2)
            plot(d_proOlamc_Near_wall,cdfd_proOlam_Near_wall,...
                'LineStyle','none','color','b','Marker','^');
            hold on
            plot(d_retOlamc_Near_wall,cdfd_retOlam_Near_wall,...
                'LineStyle','none','color','r','Marker','*');
            plot(x_values,fitted_cdf,...
                'k','LineWidth',1.5)
            set(gca,'TickLabelInterpreter','latex','FontSize',13)
            xlabel('$d_{w}/\lambda_{T}$','Interpreter','Latex','FontSize',14);
            ylabel('c.d.f.','Interpreter','Latex','FontSize',14);
            legend('Near wall Prograge','Near wall Retrograde','Gen','Interpreter','latex','FontSize',9)
%             close all
            z_prog = zeros(1,1);
            z_retro = zeros(1,1);
            pr = 1;
            rt = 1;
            for S=1:size(Obj.HRVFu,3)
                xend_near_wall_prev =0;

                while xend_near_wall_prev < max(Obj.HRVFx)
                    
                    N = 2;%2
                    r_omega = dwOlambda_T(ii)*Obj.Delx*0.5;
                    zc_near_wall = N*r_omega + Obj.HRVFz(1);
                    xc_near_wall = N*r_omega + xend_near_wall_prev;
                    xend_near_wall_prev = xc_near_wall + N*r_omega;%6
                    [indzmid] = find(Obj.HRVFz>=zc_near_wall,1);
                    [indxmid] = find(Obj.HRVFx>=xc_near_wall,1);
                    if Obj.HRVFuprime(indzmid,indxmid,S)<0
                       ii = ii + 1;
                       continue 
                    end
                    Gama=2*pi*0.5*delUwOu_tau(ii)*Obj.u_tau*r_omega;
                    ii = ii + 1;
                    [indz] = find(Obj.HRVFz>=zc_near_wall-N*r_omega &...
                        Obj.HRVFz<=zc_near_wall+N*r_omega);
                    [indx] = find(Obj.HRVFx>=xc_near_wall-N*r_omega &...
                        Obj.HRVFx<=xc_near_wall+N*r_omega);
                    [X,Z] = meshgrid(Obj.HRVFx(indx),Obj.HRVFz(indz));
                    r = sqrt((X-xc_near_wall).^2+(Z-zc_near_wall).^2);
%                     r(r>N*1.01*r_omega)=0;
                    z_prog(pr,1) = zc_near_wall;
                    pr = pr +1;
                    
                    uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
                    nan_logical_array = isnan(uazi);
                    [nan_indices] = find(nan_logical_array);
                    uazi(nan_indices) = 0;
                    
                    Theta = atan2((Z-zc_near_wall),(X-xc_near_wall));
                    nan_logical_array = isnan(Theta);
                    [nan_indices] = find(nan_logical_array);
                    Theta(nan_indices) = 0;
                    
%                     mask1 = (Theta >= 0 & Theta < pi/2 & r~=0) ;
%                     
%                     mask2 =(Theta >= pi/2 & Theta <= pi+0.1 & r~=0);
%                     
%                     mask3 =((Theta >= -pi & Theta <= -pi/2 & r~=0));
%                     
%                     mask4 = (Theta >= -pi/2 & Theta <= 0 & r~=0) ;
                    
                    
                    u_vor = uazi.*sin(Theta);
                    u_vor(uazi~=0)=u_vor(uazi~=0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                    w_vor = -uazi.*cos(Theta);
%                     w_vor(uazi>0)=w_vor(uazi>0)-Gama/(2*pi*r_omega)*(1-exp(-1));
%                     A1 = zeros(size(w_vor(mask1),1),2);
%                     A2 = zeros(size(w_vor(mask2),1),2);
%                     A3 = zeros(size(w_vor(mask3),1),2);
%                     A4 = zeros(size(w_vor(mask4),1),2);
%                     A1(:,1) = u_vor(mask1);
%                     A1(:,2) = w_vor(mask1);
%                     A2(:,1) = u_vor(mask2);
%                     A2(:,2) = w_vor(mask2);
%                     A3(:,1) = u_vor(mask3);
%                     A3(:,2) = w_vor(mask3);
%                     A4(:,1) = u_vor(mask4);
%                     A4(:,2) = w_vor(mask4);
%                     
%                     
%                     B1 = A1* Trans;
%                     B2 = A2* Trans;
%                     B3 = A3* Trans;
%                     B4 = A4* Trans;
%                     u_vor_trans = zeros(size(u_vor));
%                     w_vor_trans = zeros(size(w_vor));
%                     u_vor_trans(mask1) = B1(:,1);
%                     w_vor_trans(mask1) = B1(:,2);
%                     
%                     u_vor_trans(mask2) = B2(:,1);
%                     w_vor_trans(mask2) = B2(:,2);
%                     
%                     u_vor_trans(mask3) = B3(:,1);
%                     w_vor_trans(mask3) = B3(:,2);
%                     
%                     u_vor_trans(mask4) = B4(:,1);
%                     w_vor_trans(mask4) = B4(:,2);
                    
%                     figure
%                     quiver((X(1:2:end,1:2:end)-xc_near_wall)./r_omega,...
%                         (Z(1:2:end,1:2:end)-zc_near_wall)./r_omega,...
%                         u_vor_trans(1:2:end,1:2:end),w_vor_trans(1:2:end,1:2:end)...
%                         ,2,'color',[1.00,0.00,0.00],...
%                         'MaxHeadSize',0.2,'LineWidth',1.5);
%                     set(gca,'TickLabelInterpreter','latex','FontSize',13,...
%                         'XGrid','on','YGrid','on')
%                     xlabel('x/r$_{\omega}$','Interpreter','Latex','FontSize',14);
%                     ylabel('z/r$_{\omega}$','Interpreter','Latex','FontSize',14);
%                     axis equal
%                     xlim([-1.2 1.2])
%                     ylim([-1.2 1.2])
%                     
%                     Obj.HRVFu(indz,indx,S) = Obj.HRVFu(indz,indx,S)+ u_vor_trans;
% %                     Obj.HRVFuprime(indz,indx,S) = Obj.HRVFuprime(indz,indx,S)+ u_vor_trans;
%                     Obj.HRVFw(indz,indx,S) = Obj.HRVFw(indz,indx,S)+ w_vor_trans;
                    Obj.HRVFu(indz,indx,S) = Obj.HRVFu(indz,indx,S)+ u_vor;
%                     Obj.HRVFuprime(indz,indx,S) = Obj.HRVFuprime(indz,indx,S)+ u_vor;
                    Obj.HRVFw(indz,indx,S) = Obj.HRVFw(indz,indx,S)+ w_vor;
                end
                
                %                 Obj.HRVFu(:,:,S)=0;
                %                 Obj.HRVFw(:,:,S)=0;
                %Location of the vortex now is at the center(it can be stochastic
                %or based onthe swirling motion)
                %Number of the vortices should be changed baased on the
                %number of lambda_ci
                Matrix = Obj.Lambda_ci(:,:,S);
                
                binaryMatrix = Matrix ~= 0;
                
                % Label connected components
                [labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);  % 4-connectivity
                
                % Initialize variables to store centroids
                centroids = zeros(1,5);
                ix = 1;
                % Calculate the mass center (centroid) for each connected component
                for i = 1:numComponents
                    % Extract the indices of the current component
                    [rows, cols] = find(labeledMatrix == i);
                    values = Matrix(labeledMatrix == i);
                    sign_region_mean = sign(mean(values));
                    % Calculate the centroid
                    row_centroid = sum(rows .* values) / sum(values);
                    col_centroid = sum(cols .* values) / sum(values);
                    
                    % Round the centroids to the nearest integer
                    row_centroid = round(row_centroid);
                    col_centroid = round(col_centroid);
                    if row_centroid <= 0 || row_centroid > size(Obj.HRVFz,1)
                        continue
                    end
                    if col_centroid <= 0 || col_centroid > size(Obj.HRVFx,2)
                        continue
                    end
                    % Store the centroid
                    centroids(ix, :) = [row_centroid, col_centroid, sign_region_mean,...
                        min(cols), max(cols)];
                    ix = ix + 1;
                end
                
                
                xcs = Obj.HRVFx(centroids(:, 2));
                
                zcs = Obj.HRVFz(centroids(:, 1));
                
                rot = centroids(:, 3);
                Numbvortices  = size(centroids,1);
                %                 indxcs = randi([1,size(Obj.HRVFu,2)],Numbvortices,1);
                %                 xcs = Obj.HRVFx(indxcs);
                %                 indzcs = randi([1,size(Obj.HRVFu,1)],Numbvortices,1);
                %                 zcs = Obj.HRVFz(indzcs);
                %                 xcs = Obj.HRVFx(floor(size(Obj.HRVFu,2)/2));
                %                 zcs = Obj.HRVFz(floor(size(Obj.HRVFu,1)/2));
                %                 kappa = 0.39;
                for i=1:Numbvortices
                    N1 = 2;%This should be 2 otherwise C1, C2 has no significant effect
                    C1 = 1;
                    C2 = 1;
                    [Obj,ii,z_prog, z_retro, pr, rt] = lambdacivorX_mac(Obj,i,ii,S,xcs,zcs,rot,dwOlambda_T,...
                        delUwOu_tau, z_prog, z_retro, pr, rt, Trans, N1, C1, C2);
                    xold = xcs(i);
                    N2 = 1;
                    while  xold+N1*dwOlambda_T(ii-1)*Obj.Delx*0.5+N2*dwOlambda_T(ii)*Obj.Delx*0.5<Obj.HRVFx(centroids(i, 5))
                        [Obj,ii,z_prog, z_retro, pr, rt] = lambdacivorX_mac(Obj,i,ii,S,xcs,zcs,rot,dwOlambda_T,...
                        delUwOu_tau, z_prog, z_retro, pr, rt, Trans, N2, C1, C2);
                        xold = xold+N1*dwOlambda_T(ii-1)*Obj.Delx*0.5+N2*dwOlambda_T(ii)*Obj.Delx*0.5;
                        N1 = N2;
                    end
                    xold = xcs(i);
                    N1 = 2;
                    ii = ii +1;
                    while xold-N1*dwOlambda_T(ii-1)*Obj.Delx*0.5-N2*dwOlambda_T(ii)*Obj.Delx*0.5>Obj.HRVFx(centroids(i, 4))
                        [Obj,ii,z_prog, z_retro, pr, rt] = lambdacivorX_mac(Obj,i,ii,S,xcs,zcs,rot,dwOlambda_T,...
                        delUwOu_tau, z_prog, z_retro, pr, rt, Trans, N2, C1, C2);
                        xold = xold-N1*dwOlambda_T(ii-1)*Obj.Delx*0.5-N2*dwOlambda_T(ii)*Obj.Delx*0.5;
                        N1 = N2;
                    end
                end                
            end
%             Obj.HRVFuprime = Obj.HRVFu - mean(Obj.HRVFu,3);
            Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);
            numb_prog = zeros(size(Obj.HRVFz)-1);
            numb_retro = zeros(size(Obj.HRVFz)-1);
            for k = 1:size(Obj.HRVFz,1)-1
                numb_prog(k,1) = sum(z_prog>=Obj.HRVFz(k,1)&...
                    z_prog<=Obj.HRVFz(k+1,1),1);
                numb_retro(k,1) = sum(z_retro>=Obj.HRVFz(k,1)&...
                    z_retro<=Obj.HRVFz(k+1,1),1);
            end
            figure
            subplot(1,2,1)
            plot(numb_prog/(Obj.HRVFx(end)^2*Obj.HRVFz(end)^2*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.delta))
            set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
                'YGrid','on')
            ylabel('z/$\delta$','Interpreter','Latex','FontSize',14);
            xlabel('Number of prograde/(Area)','Interpreter','Latex','FontSize',14);
            subplot(1,2,2)
            plot(numb_retro/(Obj.HRVFx(end)^2*Obj.HRVFz(end)^2*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.delta),'r')
            set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
                'YGrid','on')
            ylabel('z/$\delta$','Interpreter','Latex','FontSize',14);
            xlabel('Number of retrograde/(Area)','Interpreter','Latex','FontSize',14);
        end
        function Newfunction(Obj, propName)
            if ~isprop(Obj, propName)
                addprop(Obj, propName);
            end
            Obj.(propName) = zeros(size(Obj.u));
            
        end
        
    end
end