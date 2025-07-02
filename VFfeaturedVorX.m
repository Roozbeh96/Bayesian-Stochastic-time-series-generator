classdef VFfeaturedVorX < dynamicprops
    properties
        name, type, u, uprime, w, x, z, delta, u_tau
        Lambda_ci, LSM, Delx, HRVFu, HRVFw, HRVFx
        HRVFz, factor, HRVFuprime, Rij, Lambda_cirms, z0, nu
        PowSpecDen_k, eta_spec, epsilon_spec, wavenumb, fs,
        D11,r11,eta_str, epsilon_str, ks, Delta_u_m, z_shear, x_shear, Lambda_T
        %dudx, dudz
        %dwdx, dwdz
    end

    methods
        function Obj = VFfeaturedVorX(name, type, u, uprime, w, x, z, delta, u_tau,...
                LSM, Delx, z0, nu, fs)

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
            Obj.z0 = z0;
            Obj.nu = nu;
            Obj.fs = fs;
            Obj.ks = z0/exp(-0.39*8.5);
        end
        function detecting_shearlayers(Obj)
            Obj.Delta_u_m = cell(size(Obj.HRVFu,2),size(Obj.HRVFu,3));
            Obj.z_shear = cell(size(Obj.HRVFu,2),size(Obj.HRVFu,3));
            Obj.x_shear = cell(size(Obj.HRVFu,2),size(Obj.HRVFu,3));
            for S = 1:size(Obj.HRVFu,3)
                for p=1:size(Obj.HRVFu,2)

                    [unique_values, first_indices] = unique(Obj.HRVFu(:,p,S), 'stable');
                    Obj.Delta_u_m{p,S} = cell(1,length(unique_values)-1);
                    Obj.z_shear{p,S} = cell(1,length(unique_values)-1);
                    Obj.x_shear{p,S} = cell(1,length(unique_values)-1);
                    if length(unique_values) == 1
                        continue
                    end
                    for umz=1:length(unique_values)-1
                        Obj.Delta_u_m{p,S}{1,umz} = unique_values(umz+1)-unique_values(umz);
                        Obj.z_shear{p,S}{1,umz} = Obj.z(first_indices(umz+1));
                        Obj.x_shear{p,S}{1,umz} = Obj.HRVFx(p);
                    end

                end
            end
        end
        function Increasing_Res_nonlinear(Obj, factor)
            if Obj.type == "Exp"

            else
                L = ceil(Obj.LSM*Obj.delta/Obj.Delx);%Number of intervals
                Numrator = floor((size(Obj.u,2)-1)/L);
                DelxGen = Obj.Delx;
                DelzGen = Obj.z(2)-Obj.z(1);
                kappa = 0.39;
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

                for S= 1:Numrator

                    ucroped = Obj.u(:, (S-1)*L+1:(S)*L+1);
                    wcroped = Obj.w(:, (S-1)*L+1:(S)*L+1);

                    for r = 1:size(Obj.HRVFz,1)
                        % Interpolate the row
                        Obj.HRVFu(r, :, S) = interp1(xorg, ucroped(r, :), xi, 'makima');
                        Obj.HRVFw(r, :, S) = interp1(xorg, wcroped(r, :), xi, 'makima');
                    end
                end
                %                 Obj.HRVFu = repmat(Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0), [1, size(Obj.HRVFu,2), size(Obj.HRVFu,3)]); %just for ASL if we want to use sto_gen uprime for mean flow
                %                 Obj.HRVFuprime = Obj.HRVFu-mean(Obj.HRVFu,3);
                Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);


            end

        end
        function Increasing_Res_sto(Obj, factor)
            if Obj.type == "Exp"

            else
                %                 if Obj.name == "GenWT7"
%                 std_uprime_shorttime = load('std_uprime_shorttime_short_GenWT7.mat').std_uprime_shorttime_GenWT7;
%                 std_wprime_shorttime = load('std_wprime_shorttime_short_GenWT7.mat').std_wprime_shorttime_GenWT7;

                std_uprime_shorttime = load('std_uprime_shorttime_Long_GenWT7.mat').std_uprime_shorttime_GenWT7;
                std_wprime_shorttime = load('std_wprime_shorttime_Long_GenWT7.mat').std_wprime_shorttime_GenWT7;

                %                     rho = -0.3;
                %                 elseif Obj.name == "GenWT10"
                %                     std_uprime_shorttime = load('std_uprime_shorttime_GenWT10.mat').std_uprime_shorttime_GenWT10;
                %                     std_wprime_shorttime = load('std_wprime_shorttime_GenWT10.mat').std_wprime_shorttime_GenWT10;
                %                     rho = -0.3;
                %                 else
                %                     std_uprime_shorttime = load('std_uprime_shorttime_GenASL.mat').std_uprime_shorttime_GenASL;
                %                     std_wprime_shorttime = load('std_wprime_shorttime_GenASL.mat').std_wprime_shorttime_GenASL;
                %                     std_uprime_shorttime = std_uprime_shorttime;
                %                     std_wprime_shorttime = std_wprime_shorttime;
                %                     rho = -0.22;
                %                 end




                %                 for i=1:size(Given_dist_NLmin_NLmax,2)
                %                     [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                %                         "ReturnList",x_test = Given_dist_NLmin_NLmax(i))
                %                     mix_coeff_Given_dist_NLmin_NLmax(i,:) = result{1};
                %                     mean_mix_Given_dist_NLmin_NLmax(i,:) = result{2};
                %                     std_mix_Given_dist_NLmin_NLmax(i,:) = result{3};
                %                 end


                rho = -0.3;
                L = ceil(Obj.LSM*Obj.delta/Obj.Delx);%Number of intervals
                Numrator = floor((size(Obj.u,2)-1)/L);
                DelxGen = Obj.Delx;
                DelzGen = Obj.z(2)-Obj.z(1);
                kappa = 0.39;
                if isnan(factor)
                    Obj.factor = floor(DelxGen/DelzGen);
                else
                    Obj.factor = factor;
                end
                %                 Obj.HRVFu=zeros(size(Obj.u,1), Obj.factor*(L) +1, Numrator);
                %                 Obj.HRVFw=zeros(size(Obj.u,1), Obj.factor*(L) +1, Numrator);
                %                 Obj.HRVFx= 0:DelxGen/Obj.factor: (L)*DelxGen;
                %                 Obj.HRVFz=Obj.z;
                xorg = 0:DelxGen:(L)*DelxGen;
                Obj.x = xorg;
                %                 xi = Obj.HRVFx;
                xi = size(Obj.HRVFu,2)*size(Obj.HRVFu,3);
                res = (Obj.HRVFx(2)-Obj.HRVFx(1));

%                 size_minmaxvec = 2.5e6;
%                 minmaxvec = zeros(size_minmaxvec,2);
%                 minind = 1;
%                 maxind = 1;
%                 start = std_uprime_shorttime(1,1)*randn(1,1)/Obj.u_tau;
%                 sample = std_uprime_shorttime(1,1)*randn(1,1)/Obj.u_tau;
%                 if start>sample
%                     minmaxvec(minind,1) = sample;
%                     minind = minind+1;
%                 else
%                     minmaxvec(maxind,2) = sample;
%                     maxind = maxind+1;
%                 end
%                 is_min_sampled = (start > sample);  % true if start > sample, otherwise false
% 
%                 for x_ = 1:2*size_minmaxvec-1
%                     if is_min_sampled
% 
% 
%                         % Run model for the min value
%                         [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_min_model.py", ...
%                             "ReturnList", x_test = minmaxvec(minind-1,1));
%                         mix_coeff_Given_min = single(result{1});
%                         mean_mix_Given_min = single(result{2});
%                         std_mix_Given_min = single(result{3});
%                         minmaxvec(maxind,2) = minmaxvec(minind-1,1);
%                         while minmaxvec(maxind,2)<=minmaxvec(minind-1,1)
%                             minmaxvec(maxind,2)=0;
%                             randomvar =rand(1,1);
%                             for j =1:size(single(result{1}),2)
%                                 minmaxvec(maxind,2) = mix_coeff_Given_min(1,j)*norminv(randomvar,mean_mix_Given_min(1,j),std_mix_Given_min(1,j)) + minmaxvec(maxind,2);
%                             end
%                         end
% 
%                         maxind = maxind + 1;
% 
%                     else
% 
% 
%                         % Run model for the max value
%                         [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_max_model.py", ...
%                             "ReturnList", x_test = minmaxvec(maxind-1,2));
%                         mix_coeff_Given_max = single(result{1});
%                         mean_mix_Given_max = single(result{2});
%                         std_mix_Given_max = single(result{3});
%                         minmaxvec(minind,1) = minmaxvec(maxind-1,2);
%                         while minmaxvec(minind,1)>=minmaxvec(maxind-1,2)
%                             minmaxvec(minind,1)=0;
%                             randomvar =rand(1,1);
%                             for j =1:size(single(result{1}),2)
%                                 minmaxvec(minind,1) = mix_coeff_Given_max(1,j)*norminv(randomvar,mean_mix_Given_max(1,j),std_mix_Given_max(1,j)) + minmaxvec(minind,1);
%                             end
%                         end
%                         minind = minind + 1;
%                     end
% 
%                     % Toggle the flag after each iteration
%                     is_min_sampled = ~is_min_sampled;
%                 end
%                 save('minmaxvec_short.mat','minmaxvec')
%                 data = load('minmaxvec_short.mat');
                data = load('minmaxvec.mat');
                minmaxvec = data.minmaxvec;
                minind = 1;
                maxind = 1;
                for r = 1:size(Obj.HRVFz,1)
                    if r ==1


                        SGuprime = zeros(size(Obj.z,1),size(Obj.HRVFu,2)*size(Obj.HRVFu,3));
                        SGwprime = zeros(size(Obj.z,1),size(Obj.HRVFu,2)*size(Obj.HRVFu,3));
                        %                             minind = 1;
                        %                             maxind = 1;
                        start = std_uprime_shorttime(r,1)*randn(1,1)/Obj.u_tau;
                        SGuprime(1,1) = start;
                        %                             sample = std_uprime_shorttime(r,1)*randn(1,1)/Obj.u_tau;

                        if start >= minmaxvec(1,1)
                            is_min_sampled = true;
                        else
                            is_min_sampled = false;
                        end

                        first_toward = true;
                        i=2;
                        while i <= xi
                            if first_toward
                                if is_min_sampled

                                    duNLmin = abs(SGuprime(r,i-1)-minmaxvec(minind,1));
                                    [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                                        "ReturnList",x_test = duNLmin);
                                    mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                                    mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                                    std_mix_Given_dist_NLmin_NLmax = single(result{3});
                                    randomvar =rand(1,1);
                                    duNN = 0;
                                    for j =1:size(single(result{1}),2)
                                        duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                                            std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                                    end


                                    if (SGuprime(r,i-1)-exp(duNN)*res)<minmaxvec(minind,1)

                                        first_toward = false;
                                        SGuprime(r,i) = minmaxvec(minind,1);
                                        is_min_sampled =~is_min_sampled;
                                        minind = minind + 1;
                                        i = i+1;

                                    else
                                        SGuprime(r,i)=SGuprime(r,i-1)-exp(duNN)*res;
                                        i = i+1;
                                    end

                                else
                                    duNLmax = abs(SGuprime(r,i-1)-minmaxvec(maxind,2));

                                    [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                                        "ReturnList",x_test = duNLmax);
                                    mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                                    mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                                    std_mix_Given_dist_NLmin_NLmax = single(result{3});
                                    randomvar =rand(1,1);
                                    duNN = 0;
                                    for j =1:size(single(result{1}),2)
                                        duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                                            std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                                    end


                                    if (SGuprime(r,i-1)+exp(duNN)*res) > minmaxvec(maxind,2)

                                        first_toward = false;
                                        SGuprime(r,i) = minmaxvec(maxind,2);
                                        is_min_sampled =~is_min_sampled;
                                        maxind = maxind + 1;
                                        i = i+1;

                                    else
                                        SGuprime(r,i)=SGuprime(r,i-1)+exp(duNN)*res;
                                        i=i+1;
                                    end
                                end
                            else
                                if is_min_sampled
                                    duNLmax = abs(SGuprime(r,i-1)-minmaxvec(maxind-1,2));
                                    duNLmin = abs(SGuprime(r,i-1)-minmaxvec(minind,1));
                                    [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                                        "ReturnList",x_test = min(duNLmin,duNLmax));
                                    mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                                    mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                                    std_mix_Given_dist_NLmin_NLmax = single(result{3});
                                    randomvar =rand(1,1);
                                    duNN = 0;
                                    for j =1:size(single(result{1}),2)
                                        duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                                            std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                                    end


                                    if (SGuprime(r,i-1)-exp(duNN)*res)<minmaxvec(minind,1)


                                        SGuprime(r,i) = minmaxvec(minind,1);
                                        is_min_sampled =~is_min_sampled;
                                        minind = minind + 1;
                                        i = i+1;

                                    else
                                        SGuprime(r,i)=SGuprime(r,i-1)-exp(duNN)*res;
                                        i = i+1;
                                    end

                                else
                                    duNLmax = abs(SGuprime(r,i-1)-minmaxvec(maxind,2));
                                    duNLmin = abs(SGuprime(r,i-1)-minmaxvec(minind-1,1));
                                    [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                                        "ReturnList",x_test = min(duNLmax,duNLmin));
                                    mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                                    mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                                    std_mix_Given_dist_NLmin_NLmax = single(result{3});
                                    randomvar =rand(1,1);
                                    duNN = 0;
                                    for j =1:size(single(result{1}),2)
                                        duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                                            std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                                    end


                                    if (SGuprime(r,i-1)+exp(duNN)*res) > minmaxvec(maxind,2)


                                        SGuprime(r,i) = minmaxvec(maxind,2);
                                        is_min_sampled =~is_min_sampled;
                                        maxind = maxind + 1;
                                        i = i+1;

                                    else
                                        SGuprime(r,i)=SGuprime(r,i-1)+exp(duNN)*res;
                                        i=i+1;
                                    end
                                end

                            end
                        end

                        sigma1 = 1*std_uprime_shorttime(r,1);
                        sigma2 = 1*std_wprime_shorttime(r,1);

                        SGuprime(r,:) = SGuprime(r,:) * Obj.u_tau;
                        SGuprime(r,:) = SGuprime(r,:) - mean(SGuprime(r,:),2);

                        SGwprime(r,:) = rho* sigma2/sigma1 * SGuprime(r,:) + sqrt(1 - rho^2) * normrnd(0, sigma2, [1, size(xi,2)]);
                        SGwprime(r,:) = SGwprime(r,:) - mean(SGwprime(r,:),2);

                    else
                        sigma1 = 1*std_uprime_shorttime(r,1);
                        sigma2 = 1*std_wprime_shorttime(r,1);

                        SGuprime(r,:) = 0.99 * SGuprime(r-1,:) + sqrt(1 - 0.99^2) * normrnd(0, sigma1, [1, size(xi,2)]);
                        SGuprime(r,:) = SGuprime(r,:) - mean(SGuprime(r,:),2);

                        SGwprime(r,:) = rho* sigma2/sigma1 * SGuprime(r,:) + sqrt(1 - rho^2) * normrnd(0, sigma2, [1, size(xi,2)]);
                        SGwprime(r,:) = SGwprime(r,:) - mean(SGwprime(r,:),2);
                    end

                end
                for S= 1:Numrator              
                    for r = 1:size(Obj.HRVFz,1)
                        % Interpolate the row
                        Obj.HRVFu(r, :, S) = Obj.HRVFu(r, :, S) + SGuprime(r,(S-1)*(Obj.factor*(L) +1)+1:(S)*(Obj.factor*(L) +1));
                        Obj.HRVFw(r, :, S) = Obj.HRVFw(r, :, S) + SGwprime(r,(S-1)*(Obj.factor*(L) +1)+1:(S)*(Obj.factor*(L) +1));
                    end
                end

                Obj.HRVFuprime = Obj.HRVFu-mean(Obj.HRVFu,3);
                %                 for S= 1:Numrator
                %
                %                     for r = 1:size(Obj.HRVFz,1)
                %                         if r ==1
                %
                %
                %                             SGuprime = zeros(size(Obj.z,1),size(Obj.HRVFx,2));
                % %                             minind = 1;
                % %                             maxind = 1;
                %                             start = std_uprime_shorttime(r,1)*randn(1,1)/Obj.u_tau;
                %                             SGuprime(1,1) = start;
                % %                             sample = std_uprime_shorttime(r,1)*randn(1,1)/Obj.u_tau;
                %
                %                             if minind<=maxind
                %                                 is_min_sampled = true;
                %                             else
                %                                 is_min_sampled = false;
                %                             end
                %
                %                             first_toward = true;
                %                             i=2;
                %                             while i <= size(xi,2)
                %                                 if first_toward
                %                                     if is_min_sampled
                %
                %                                         duNLmin = abs(SGuprime(r,i-1)-minmaxvec(minind,1));
                %                                         [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                %                                             "ReturnList",x_test = duNLmin);
                %                                         mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                %                                         mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                %                                         std_mix_Given_dist_NLmin_NLmax = single(result{3});
                %                                         randomvar =rand(1,1);
                %                                         duNN = 0;
                %                                         for j =1:size(single(result{1}),2)
                %                                             duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                %                                                 std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                %                                         end
                %
                %
                %                                         if (SGuprime(r,i-1)-exp(duNN)*res)<minmaxvec(minind,1)
                %
                %                                             first_toward = false;
                %                                             SGuprime(r,i) = minmaxvec(minind,1);
                %                                             is_min_sampled =~is_min_sampled;
                %                                             minind = minind + 1;
                %                                             i = i+1;
                %
                %                                         else
                %                                             SGuprime(r,i)=SGuprime(r,i-1)-exp(duNN)*res;
                %                                             i = i+1;
                %                                         end
                %
                %                                     else
                %                                         duNLmax = abs(SGuprime(r,i-1)-minmaxvec(maxind,2));
                %
                %                                         [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                %                                             "ReturnList",x_test = duNLmax);
                %                                         mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                %                                         mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                %                                         std_mix_Given_dist_NLmin_NLmax = single(result{3});
                %                                         randomvar =rand(1,1);
                %                                         duNN = 0;
                %                                         for j =1:size(single(result{1}),2)
                %                                             duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                %                                                 std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                %                                         end
                %
                %
                %                                         if (SGuprime(r,i-1)+exp(duNN)*res) > minmaxvec(maxind,2)
                %
                %                                             first_toward = false;
                %                                             SGuprime(r,i) = minmaxvec(maxind,2);
                %                                             is_min_sampled =~is_min_sampled;
                %                                             maxind = maxind + 1;
                %                                             i = i+1;
                %
                %                                         else
                %                                             SGuprime(r,i)=SGuprime(r,i-1)+exp(duNN)*res;
                %                                             i=i+1;
                %                                         end
                %                                     end
                %                                 else
                %                                     if is_min_sampled
                %                                         duNLmax = abs(SGuprime(r,i-1)-minmaxvec(maxind-1,2));
                %                                         duNLmin = abs(SGuprime(r,i-1)-minmaxvec(minind,1));
                %                                         [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                %                                             "ReturnList",x_test = min(duNLmin,duNLmax));
                %                                         mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                %                                         mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                %                                         std_mix_Given_dist_NLmin_NLmax = single(result{3});
                %                                         randomvar =rand(1,1);
                %                                         duNN = 0;
                %                                         for j =1:size(single(result{1}),2)
                %                                             duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                %                                                 std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                %                                         end
                %
                %
                %                                         if (SGuprime(r,i-1)-exp(duNN)*res)<minmaxvec(minind,1)
                %
                %
                %                                             SGuprime(r,i) = minmaxvec(minind,1);
                %                                             is_min_sampled =~is_min_sampled;
                %                                             minind = minind + 1;
                %                                             i = i+1;
                %
                %                                         else
                %                                             SGuprime(r,i)=SGuprime(r,i-1)-exp(duNN)*res;
                %                                             i = i+1;
                %                                         end
                %
                %                                     else
                %                                         duNLmax = abs(SGuprime(r,i-1)-minmaxvec(maxind,2));
                %                                         duNLmin = abs(SGuprime(r,i-1)-minmaxvec(minind-1,1));
                %                                         [result] = pyrunfile("G:\My Drive\Research\VFfeaturedVorX\Given_dist_NLmin_NLmax_model.py",...
                %                                             "ReturnList",x_test = min(duNLmax,duNLmin));
                %                                         mix_coeff_Given_dist_NLmin_NLmax = single(result{1});
                %                                         mean_mix_Given_dist_NLmin_NLmax = single(result{2});
                %                                         std_mix_Given_dist_NLmin_NLmax = single(result{3});
                %                                         randomvar =rand(1,1);
                %                                         duNN = 0;
                %                                         for j =1:size(single(result{1}),2)
                %                                             duNN = mix_coeff_Given_dist_NLmin_NLmax(1,j)*norminv(randomvar,mean_mix_Given_dist_NLmin_NLmax(1,j),...
                %                                                 std_mix_Given_dist_NLmin_NLmax(1,j)) + duNN;
                %                                         end
                %
                %
                %                                         if (SGuprime(r,i-1)+exp(duNN)*res) > minmaxvec(maxind,2)
                %
                %
                %                                             SGuprime(r,i) = minmaxvec(maxind,2);
                %                                             is_min_sampled =~is_min_sampled;
                %                                             maxind = maxind + 1;
                %                                             i = i+1;
                %
                %                                         else
                %                                             SGuprime(r,i)=SGuprime(r,i-1)+exp(duNN)*res;
                %                                             i=i+1;
                %                                         end
                %                                     end
                %
                %                                 end
                %                             end
                %
                %                             sigma1 = 1*std_uprime_shorttime(r,1);
                %                             sigma2 = 1*std_wprime_shorttime(r,1);
                %
                %                             SGuprime(r,:) = SGuprime(r,:) * Obj.u_tau;
                %
                %                             SGwprime = rho* sigma2/sigma1 * SGuprime(r,:) + sqrt(1 - rho^2) * normrnd(0, sigma2, [1, size(xi,2)]);
                %
                %
                %
                %                             Obj.HRVFu(r, :, S) = Obj.HRVFu(r, :, S) + SGuprime(r,:) ;
                %                             Obj.HRVFw(r, :, S) = Obj.HRVFw(r, :, S) + SGwprime;
                %                         else
                %                             sigma1 = 1*std_uprime_shorttime(r,1);
                %                             sigma2 = 1*std_wprime_shorttime(r,1);
                %
                %                             SGuprime(r,:) = 0.99 * SGuprime(r-1,:) + sqrt(1 - 0.99^2) * normrnd(0, sigma1, [1, size(xi,2)]);
                %                             SGwprime = rho* sigma2/sigma1 * SGuprime(r,:) + sqrt(1 - rho^2) * normrnd(0, sigma2, [1, size(xi,2)]);
                %
                %                             Obj.HRVFu(r, :, S) = Obj.HRVFu(r, :, S)+ SGuprime(r,:);
                %                             Obj.HRVFw(r, :, S) = Obj.HRVFw(r, :, S)+ SGwprime;
                %                         end
                %
                %                     end
            end
            %                 Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);
            %                   Obj.HRVFw = Obj.HRVFw-mean(Obj.HRVFw,3);
        end
        function Lambdaci(Obj)
            if Obj.type == "Exp"
                Obj.Lambda_ci = zeros(size(Obj.u));
                DelxExp = Obj.x(2)-Obj.x(1);
                DelzExp = Obj.z(2)-Obj.z(1);
            else
                %             if isempty(Obj.Lambda_ci)
                Obj.Lambda_ci = zeros(size(Obj.HRVFu));
                % %                                     Obj.omega = zeros(size(Obj.HRVFu));
                %             end
                DelxGen = Obj.HRVFx(2)-Obj.HRVFx(1);
                DelzGen = Obj.z(2)-Obj.z(1);
            end
            if Obj.type == "Exp"
                for S= 1:size(Obj.Lambda_ci, 3)

                    [dudx, dudz] = gradient(Obj.u(:,:, S),DelxExp,DelzExp);
                    [dwdx, dwdz] = gradient(Obj.w(:,:, S),DelxExp,DelzExp);
                    %                     Obj.omega(:,:,S) = dwdx - dudz;
                    omega = dwdx - dudz;
                    Mat = zeros(size(dudx,1),size(dudx,2),2,2);
                    Mat(:,:,1,1) = dudx;
                    Mat(:,:,1,2) = dudz;
                    Mat(:,:,2,1) = dwdx;
                    Mat(:,:,2,2) = dwdz;
                    for r = 1:size(Obj.Lambda_ci,1)

                        for c = 1:size(Obj.Lambda_ci,2)
                            temp = reshape(Mat(r,c,:,:),[2,2]);
                            %                             Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                            %                                 *sign(Obj.omega(r,c,S));
                            Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                                *sign(omega(r,c));

                        end
                    end
                end
                Obj.Lambda_cirms = (mean(mean(Obj.Lambda_ci.^2,3),2)).^0.5;
%                 for r = 1:size(Obj.Lambda_ci, 1)
%                     Obj.Lambda_ci(r, :, :) = Obj.Lambda_ci(r, :, :) .* (abs(Obj.Lambda_ci(r, :, :)) >= 1.5*Obj.Lambda_cirms(r));
%                 end
                %                 Obj.Lambda_ci(find(abs(Obj.Lambda_ci)<1.5*Obj.Lambda_cirms))=0;
            else
                for S= 1:size(Obj.Lambda_ci,3)

                    [dudx, dudz] = gradient(Obj.HRVFu(:,:, S),DelxGen,DelzGen);
                    [dwdx, dwdz] = gradient(Obj.HRVFw(:,:, S),DelxGen,DelzGen);
                    %                     Obj.omega(:,:,S) = dwdx - dudz;
                    omega = dwdx - dudz;
                    Mat = zeros(size(dudx,1),size(dudx,2),2,2);
                    Mat(:,:,1,1) = dudx;
                    Mat(:,:,1,2) = dudz;
                    Mat(:,:,2,1) = dwdx;
                    Mat(:,:,2,2) = dwdz;
                    for r = 1:size(Obj.Lambda_ci,1)

                        for c = 1:size(Obj.Lambda_ci,2)
                            temp = reshape(Mat(r,c,:,:),[2,2]);
                            %                             Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                            %                                 *sign(Obj.omega(r,c,S));
                            Obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                                *sign(omega(r,c));

                        end
                    end


                end
                Obj.Lambda_cirms = (mean(mean(Obj.Lambda_ci.^2,3),2)).^0.5;

%                 for r = 1:size(Obj.Lambda_ci, 1)
%                     Obj.Lambda_ci(r, :, :) = Obj.Lambda_ci(r, :, :) .* (abs(Obj.Lambda_ci(r, :, :)) >= 1.0*Obj.Lambda_cirms(r));
%                 end
%                                 Obj.Lambda_ci(find(abs(Obj.Lambda_ci)<Obj.Lambda_cirms))=0;
            end

        end
        function Gauss_kernel(Obj,kern_xsize, kern_ysize, Sigma)
            kappa = 0.39;
            mu = [0 0];

            %             Sigma = [Obj.Delx*0.0005 0;...
            %                 0 Obj.Delx*0.0005];
            %             Sigma = [Obj.Delx*1 Obj.Delx*1*sqrt(tand(12.9));...
            %                 Obj.Delx*1*sqrt(tand(12.9)) tand(13)*Obj.Delx*1];
            %             Sigma = [Obj.Delx/100 sqrt(Obj.Delx/100*Obj.Delx/100*tand(0));...
            %                 sqrt(Obj.Delx/100*Obj.Delx/100*tand(0)) Obj.Delx/100];
            %             [X,Y] = meshgrid(-kern_xsize/2*Obj.Delx:Obj.HRVFx(2)-Obj.HRVFx(1):kern_xsize/2*Obj.Delx,...
            %                 -kern_ysize/2*Obj.Delx:Obj.HRVFz(2)-Obj.HRVFz(1):kern_ysize/2*Obj.Delx);
            [X,Y] = meshgrid(linspace(-kern_xsize/2 * Obj.Delx, kern_xsize/2 * Obj.Delx, ceil(kern_xsize*Obj.Delx/(Obj.HRVFx(2)-Obj.HRVFx(1)))+1),...
                linspace(-kern_ysize/2 * Obj.Delx, kern_ysize/2 * Obj.Delx, ceil(kern_ysize*Obj.Delx/(Obj.HRVFz(2)-Obj.HRVFz(1)))+1));
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
            for S =1:size(Obj.HRVFu,3)

                Obj.HRVFu(:,:,S) = conv2(Obj.HRVFu(:,:,S), kernel, 'same')./normalization_matrix;
                Obj.HRVFw(:,:,S) = conv2(Obj.HRVFw(:,:,S), kernel, 'same')./normalization_matrix;
                %                 Obj.HRVFuprime(:,:,S) = conv2(Obj.HRVFuprime(:,:,S), kernel, 'same')./normalization_matrix;
            end
            Obj.HRVFuprime = Obj.HRVFu-mean(Obj.HRVFu,3);
            %             Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);

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
        function VortXmodeling1(Obj)
            if Obj.name == "GenWT7" || Obj.name == "GenASL"
                WT7VorX = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_07ms_vortex_properties.mat');


                m = 2;

                Near_wall_dU_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*Obj.ks);
                Near_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*Obj.ks);


                Far_wall_dU_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro > m*Obj.ks);
                Far_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro > m*Obj.ks);

                R1 = corrcoef(0.5*Near_wall_d_pro_WT7,-Near_wall_dU_pro_WT7);
                rho_near = abs(R1(1,2)); % Target (Spearman) correlation %

                R2 = corrcoef(0.5*Far_wall_d_pro_WT7,-Far_wall_dU_pro_WT7);
                rho_Far = abs(R2(1,2)); % Target (Spearman) correlation %
            elseif Obj.name == "GenWT10"
                m = 2;
                WT10VorX = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_10ms_vortex_properties.mat');

                Near_wall_dU_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*Obj.ks);
                Near_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*Obj.ks);

                Far_wall_dU_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro > m*Obj.ks);
                Far_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro > m*Obj.ks);


                R1 = corrcoef(0.5*Near_wall_d_pro_WT10,-Near_wall_dU_pro_WT10);
                rho_near = abs(R1(1,2)); % Target (Spearman) correlation %

                R2 = corrcoef(0.5*Far_wall_d_pro_WT10,-Far_wall_dU_pro_WT10);
                rho_Far = abs(R2(1,2)); % Target (Spearman) correlation %
            end

            delUwOu_tau_near = load('delUwOu_tau_near.mat').delUwOu_tau_near;
            delUwOu_tau_Far = load('delUwOu_tau_Far.mat').delUwOu_tau_Far;
            rwOlambda_T_near = load('rwOlambda_T_near.mat').rwOlambda_T_near;
            rwOlambda_T_Far = load('rwOlambda_T_Far.mat').rwOlambda_T_Far;
            %             if Obj.name == "GenASL"
            %                 delUwOu_tau_near = 9*delUwOu_tau_near;
            %                 delUwOu_tau_Far = 9*delUwOu_tau_Far;
            %             end
            ii = 1;
            iii = 1;
            z_prog = zeros(1,1);
            z_retro = zeros(1,1);
            z_prog_near = zeros(1,1);
            pr = 1;
            rt = 1;
            pr_near = 1;
            rho_Near_wall = load('samples_Near.mat').samples_Near;
            rho_Far_wall = load('samples_Far.mat').samples_Far;
            xi = 1;
            xii = 1;
            Dist_close_wall_VorX = zeros(1);
            ind_close_wall_VorX = 1;
            Lambda_ci_temp = Obj.Lambda_ci;
            for r = 1:size(Obj.Lambda_ci, 1)
                    Lambda_ci_temp(r, :, :) = Obj.Lambda_ci(r, :, :) .* (abs(Obj.Lambda_ci(r, :, :)) >= 1.0*Obj.Lambda_cirms(r));
            end

            for S=1:size(Obj.HRVFu,3)
                xend_near_wall_prev = 0;
                temp_ind_close_wall_VorX = 1;
                while xend_near_wall_prev < max(Obj.HRVFx)% Attached wall vortices

                    D1 = 1;
                    D2 = 1;
                    N = 1;
                    fac = 1.0;% should be more than 1.0
                    if Obj.name ~= "GenASL"
                        r_omega = rwOlambda_T_near(ii)*Obj.Delx; %choose from near wall distribution
%                         r_omega = rwOlambda_T_near(ii)*mean(Obj.Lambda_T); %choose from near wall distribution
                        Trans = [1, rho_Near_wall(xi); 0, sqrt(1 - rho_Near_wall(xi)^2)];
                        Gama=2*pi*delUwOu_tau_near(ii)*Obj.u_tau*r_omega/(1-exp(-1)); 
                        xi= xi+1;
                        ii = ii+1;
                    else
                        r_omega = rwOlambda_T_Far(iii)*Obj.Delx;
%                         r_omega = rwOlambda_T_Far(iii)*Obj.Lambda_T;
                        Trans = [1, rho_Far_wall(xii); 0, sqrt(1 - rho_Far_wall(xii)^2)];
                        Gama=2*pi*delUwOu_tau_Far(iii)*Obj.u_tau*r_omega/(1-exp(-1)); 
                        xii = xii+1;
                        iii = iii+1;
                    end
                    zc_near_wall = N*r_omega + Obj.HRVFz(1);
                    xc_near_wall = N*r_omega + xend_near_wall_prev;
                    xend_near_wall_prev = xc_near_wall + fac*N*r_omega;
                    [indzmid] = find(Obj.HRVFz>=zc_near_wall,1);
                    [indxmid] = find(Obj.HRVFx>=xc_near_wall,1);


                    %                     logic = randi([0, 1]);

                    if Obj.HRVFuprime(indzmid,indxmid,S)<0
%                         ii = ii + 1;
%                         xi= xi+1;
%                         xii = xii+1;
                        continue
                    end

                    if randi([0, 1])==0
%                         ii = ii + 1;
%                         xi= xi+1;
%                         xii = xii+1;
                        continue
                    end

                    if temp_ind_close_wall_VorX > 1
                        Dist_close_wall_VorX(ind_close_wall_VorX-1) = xc_near_wall-prev_pos_xc_near_wall_VorX;
                    end
                    prev_pos_xc_near_wall_VorX = xc_near_wall;
                    temp_ind_close_wall_VorX = temp_ind_close_wall_VorX +1;
                    ind_close_wall_VorX = ind_close_wall_VorX +1;

                    %                     if zc_near_wall<= 2*Obj.ks

                    %                     Trans = [1, rho_Near_wall(xi); 0, sqrt(1 - rho_Near_wall(xi)^2)];
                    %                     xi= xi+1;
                    %                     else

                    %                         Trans = [1, rho_Far_wall(xii); 0, sqrt(1 - rho_Far_wall(xii)^2)];
                    %                         xii = xii+1;
                    %                     end

%                     Gama=2*pi*delUwOu_tau_near(ii)*Obj.u_tau*r_omega/(1-exp(-1));  %choose from near wall distribution
%                     ii = ii + 1;
                    [indz] = find(Obj.HRVFz>=zc_near_wall-N*r_omega &...
                        Obj.HRVFz<=zc_near_wall+N*r_omega);
                    [indx] = find(Obj.HRVFx>=xc_near_wall-N*r_omega &...
                        Obj.HRVFx<=xc_near_wall+N*r_omega);
                    [X,Z] = meshgrid(Obj.HRVFx(indx),Obj.HRVFz(indz));
                    r = sqrt((X-xc_near_wall).^2+(Z-zc_near_wall).^2);
                    %                     r(r>N*1.01*r_omega)=0;
                    z_prog_near(pr_near,1) = zc_near_wall;

                    pr_near = pr_near +1;

                    uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
                    nan_logical_array = isnan(uazi);
                    [nan_indices] = find(nan_logical_array);
                    uazi(nan_indices) = 0;

                    Theta = atan2((Z-zc_near_wall),(X-xc_near_wall));
                    nan_logical_array = isnan(Theta);
                    [nan_indices] = find(nan_logical_array);
                    Theta(nan_indices) = 0;

                    mask1 = (Theta >= 0 & Theta < pi/2 & r~=0) ;

                    mask2 =(Theta >= pi/2 & Theta <= pi+0.1 & r~=0);

                    mask3 =((Theta >= -pi & Theta <= -pi/2 & r~=0));

                    mask4 = (Theta >= -pi/2 & Theta <= 0 & r~=0) ;


                    u_vor = uazi.*sin(Theta);
                    u_vor(uazi~=0)=u_vor(uazi~=0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                    w_vor = -uazi.*cos(Theta);
                    %                     w_vor(uazi>0)=w_vor(uazi>0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                    A1 = zeros(size(w_vor(mask1),1),2);
                    A2 = zeros(size(w_vor(mask2),1),2);
                    A3 = zeros(size(w_vor(mask3),1),2);
                    A4 = zeros(size(w_vor(mask4),1),2);
                    A1(:,1) = u_vor(mask1);
                    A1(:,2) = w_vor(mask1);
                    A2(:,1) = u_vor(mask2);
                    A2(:,2) = w_vor(mask2);
                    A3(:,1) = u_vor(mask3);
                    A3(:,2) = w_vor(mask3);
                    A4(:,1) = u_vor(mask4);
                    A4(:,2) = w_vor(mask4);


                    B1 = A1* Trans;
                    B2 = A2* Trans;
                    B3 = A3* Trans;
                    B4 = A4* Trans;
                    u_vor_trans = zeros(size(u_vor));
                    w_vor_trans = zeros(size(w_vor));
                    u_vor_trans(mask1) = B1(:,1);
                    w_vor_trans(mask1) = B1(:,2);

                    u_vor_trans(mask2) = B2(:,1);
                    w_vor_trans(mask2) = B2(:,2);

                    u_vor_trans(mask3) = B3(:,1);
                    w_vor_trans(mask3) = B3(:,2);

                    u_vor_trans(mask4) = B4(:,1);
                    w_vor_trans(mask4) = B4(:,2);

                    %                     figure
                    %                     plot(B1(:,1)/(delUwOu_tau_near(ii)*Obj.u_tau),B1(:,2)/(delUwOu_tau_near(ii)*Obj.u_tau),'r.')
                    %                     hold on
                    %                     plot(B2(:,1)/(delUwOu_tau_near(ii)*Obj.u_tau),B2(:,2)/(delUwOu_tau_near(ii)*Obj.u_tau),'r.')
                    %                     plot(B3(:,1)/(delUwOu_tau_near(ii)*Obj.u_tau),B3(:,2)/(delUwOu_tau_near(ii)*Obj.u_tau),'r.')
                    %                     plot(B4(:,1)/(delUwOu_tau_near(ii)*Obj.u_tau),B4(:,2)/(delUwOu_tau_near(ii)*Obj.u_tau),'r.')
                    %                     set(gca,'TickLabelInterpreter','latex','FontSize',13,...
                    %                         'XGrid','on','YGrid','on')
                    %                     xlabel('u/$u_{\omega}$','Interpreter','Latex','FontSize',14);
                    %                     ylabel('w/$u_{\omega}$','Interpreter','Latex','FontSize',14);
                    %                     axis equal
                    %                     xlim([-1.1 1.1])
                    %                     ylim([-1.1 1.1])
                    %
                    %                     figure
                    %                     quiver((X(1:2:end,1:2:end)-xc_near_wall)./r_omega,...
                    %                         (Z(1:2:end,1:2:end)-zc_near_wall)./r_omega,...
                    %                         u_vor_trans(1:2:end,1:2:end),w_vor_trans(1:2:end,1:2:end)...
                    %                         ,1,'color',[0.86,0.09,0.24],...
                    %                         'MaxHeadSize',0.2,'LineWidth',1.5);
                    %                     set(gca,'TickLabelInterpreter','latex','FontSize',13,...
                    %                         'XGrid','on','YGrid','on')
                    %                     xlabel('x/r$_{\omega}$','Interpreter','Latex','FontSize',14);
                    %                     ylabel('z/r$_{\omega}$','Interpreter','Latex','FontSize',14);
                    %                     axis equal
                    %                     xlim([-1.2 1.2])
                    %                     ylim([-1.2 1.2])
                    if isempty(u_vor)
                        continue
                    end

                    Obj.HRVFu(indz,indx,S) = Obj.HRVFu(indz,indx,S)+ D1*u_vor_trans;
                    %                     Obj.HRVFuprime(indz,indx,S) = Obj.HRVFuprime(indz,indx,S)+ u_vor_trans;
                    Obj.HRVFw(indz,indx,S) = Obj.HRVFw(indz,indx,S)+ D2*w_vor_trans;


                    %                     Obj.HRVFu(indz,indx,S) = Obj.HRVFu(indz,indx,S)+ D1*u_vor;
                    % %                     Obj.HRVFuprime(indz,indx,S) = Obj.HRVFuprime(indz,indx,S)+ u_vor;
                    %                     Obj.HRVFw(indz,indx,S) = Obj.HRVFw(indz,indx,S)+ D2*w_vor;
                end
                % Lambda_ci vortices
                %Location of the vortex now is at the center(it can be stochastic
                %or based onthe swirling motion)
                %Number of the vortices should be changed based on the
                %number of lambda_ci
                Matrix = Lambda_ci_temp(:,:,S);

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
                    if zcs(i)<= 2*Obj.ks
                        N1 = 2;
                    else
                        N1 = 2;
                    end
                    %2This should be 2 otherwise E1, E2 has no significant effect
                    E1 = 1;
                    E2 = 1;
                    [Obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                        lambdacivorX(Obj,i,ii,iii,S,xcs,zcs,rot,rwOlambda_T_near,rwOlambda_T_Far,...
                        delUwOu_tau_near,delUwOu_tau_Far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                        rho_Far_wall(xii), xii, rho_near, rho_Far);%Primary vortex
                    xold = xcs(i);
                    if zcs(i)<= 2*Obj.ks
                        N2 = 2;
                    else
                        N2 = 2;
                    end
                    %filling periphery of the vortex(Secondary vortices)
                    if zcs(i)<= 2*Obj.ks % close to the wall
                        while  xold+N1*rwOlambda_T_near(ii-1)*Obj.Delx+N2*rwOlambda_T_near(ii)*Obj.Delx<Obj.HRVFx(centroids(i, 5))
                            [Obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                                lambdacivorX(Obj,i,ii,iii,S,xcs,zcs,rot,rwOlambda_T_near,rwOlambda_T_Far,...
                                delUwOu_tau_near,delUwOu_tau_Far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                                rho_Far_wall(xii), xii, rho_near, rho_Far);
                            xold = xold+N1*rwOlambda_T_near(ii-1)*Obj.Delx+N2*rwOlambda_T_near(ii)*Obj.Delx;
                            N1 = N2;
                        end
                        xold = xcs(i);
                        N1 = 2;%1
                        ii = ii +1;
                        while xold-N1*rwOlambda_T_near(ii-1)*Obj.Delx-N2*rwOlambda_T_near(ii)*Obj.Delx>Obj.HRVFx(centroids(i, 4))
                            [Obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                                lambdacivorX(Obj,i,ii,iii,S,xcs,zcs,rot,rwOlambda_T_near,rwOlambda_T_Far,...
                                delUwOu_tau_near,delUwOu_tau_Far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                                rho_Far_wall(xii), xii,rho_near, rho_Far);
                            xold = xold-N1*rwOlambda_T_near(ii-1)*Obj.Delx-N2*rwOlambda_T_near(ii)*Obj.Delx;
                            N1 = N2;
                        end
                    else %Far from the wall
                        while  xold+N1*rwOlambda_T_Far(iii-1)*Obj.Delx+N2*rwOlambda_T_Far(iii)*Obj.Delx<Obj.HRVFx(centroids(i, 5))
                            [Obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                                lambdacivorX(Obj,i,ii,iii,S,xcs,zcs,rot,rwOlambda_T_near,rwOlambda_T_Far,...
                                delUwOu_tau_near,delUwOu_tau_Far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                                rho_Far_wall(xii), xii,rho_near, rho_Far);
                            xold = xold+N1*rwOlambda_T_Far(iii-1)*Obj.Delx+N2*rwOlambda_T_Far(iii)*Obj.Delx;
                            N1 = N2;
                        end
                        xold = xcs(i);
                        N1 = 2;%1
                        iii = iii +1;
                        while xold-N1*rwOlambda_T_Far(iii-1)*Obj.Delx-N2*rwOlambda_T_Far(iii)*Obj.Delx>Obj.HRVFx(centroids(i, 4))
                            [Obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                                lambdacivorX(Obj,i,ii,iii,S,xcs,zcs,rot,rwOlambda_T_near,rwOlambda_T_Far,...
                                delUwOu_tau_near,delUwOu_tau_Far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                                rho_Far_wall(xii), xii,rho_near, rho_Far);
                            xold = xold-N1*rwOlambda_T_Far(iii-1)*Obj.Delx-N2*rwOlambda_T_Far(iii)*Obj.Delx;
                            N1 = N2;
                        end
                    end
                end
            end
            Obj.HRVFuprime = Obj.HRVFu - mean(Obj.HRVFu,3);
            Obj.HRVFw = Obj.HRVFw - mean(Obj.HRVFw,3);
            %             Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);
            numb_prog = zeros(size(Obj.HRVFz)-1);
            numb_prog_near = zeros(size(Obj.HRVFz)-1);
            numb_retro = zeros(size(Obj.HRVFz)-1);
            for k = 1:size(Obj.HRVFz,1)-1
                numb_prog(k,1) = sum(z_prog>=Obj.HRVFz(k,1)&...
                    z_prog<=Obj.HRVFz(k+1,1),1);
                numb_prog_near(k,1) = sum(z_prog_near>=Obj.HRVFz(k,1)&...
                    z_prog_near<=Obj.HRVFz(k+1,1),1);
                numb_retro(k,1) = sum(z_retro>=Obj.HRVFz(k,1)&...
                    z_retro<=Obj.HRVFz(k+1,1),1);
            end
            figure
            set(gcf,'Position',[941,327,806,463])
            axes('Position',[0.064516129032258,0.138228941684665,0.26302729528536,0.850971922246223])
            plot(numb_prog/(Obj.HRVFx(end)*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.z0), 'r', 'LineWidth', 2)
            set(gca,'TickLabelInterpreter','latex','FontSize',15,'XGrid','on',...
                'YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale','log')
            ylabel('z/$\mathrm{z}_{0}$','Interpreter','Latex','FontSize',15);
            xlabel('\# Prograde vortices/$\delta$','Interpreter','Latex','FontSize',15);

            axes('Position',[0.398263027295286,0.138228941684665,0.263027295285361,0.850971922246223])
            plot(numb_retro/(Obj.HRVFx(end)*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.z0),'b', 'LineWidth', 2)
            set(gca,'TickLabelInterpreter','latex','FontSize',15,'XGrid','on',...
                'YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale','log')
            ylabel('z/$\mathrm{z}_{0}$','Interpreter','Latex','FontSize',15);
            xlabel('\# Retrograde vortices/$\delta$','Interpreter','Latex','FontSize',15);


            axes('Position',[0.730769230769233,0.138228941684665,0.263027295285362,0.850971922246223])
            plot(numb_prog_near/(Obj.HRVFx(end)*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.z0), 'k', 'LineWidth', 2)
            set(gca,'TickLabelInterpreter','latex','FontSize',15,'XGrid','on',...
                'YGrid','on','XMinorGrid','on','YMinorGrid','on','YScale','log')
            ylabel('z/$\mathrm{z}_{0}$','Interpreter','Latex','FontSize',15);
            xlabel('\# Prograde near-wall vortices/$\delta$','Interpreter','Latex','FontSize',10);
            disp(['The mean value of vortex distance normalized with ks is: ', num2str(mean(Dist_close_wall_VorX,2)/(Obj.ks))]);
            disp(['The std value of vortex distance normalized with ks is: ', num2str(std(Dist_close_wall_VorX,0,2)/(Obj.ks))]);
            
        end
        function VortXmodeling2(Obj)
            if Obj.name == "GenWT7" || Obj.name == "GenASL"
                WT7VorX = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_07ms_vortex_properties.mat');


                % Find Near-wall vortices
                m = 2;

                Near_wall_dU_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*Obj.ks);
                Near_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro>=0 & WT7VorX.z_pro <= m*Obj.ks);


                Far_wall_dU_pro_WT7 = WT7VorX.dU_pro (WT7VorX.z_pro > m*Obj.ks);
                Far_wall_d_pro_WT7 = WT7VorX.d_pro (WT7VorX.z_pro > m*Obj.ks);

                R1 = corrcoef(0.5*Near_wall_d_pro_WT7,-Near_wall_dU_pro_WT7);
                rho_near = abs(R1(1,2)); % Target (Spearman) correlation %

                R2 = corrcoef(0.5*Far_wall_d_pro_WT7,-Far_wall_dU_pro_WT7);
                rho_Far = abs(R2(1,2)); % Target (Spearman) correlation %
            elseif Obj.name == "GenWT10"
                WT10VorX = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_10ms_vortex_properties.mat');
                m = 2;
                Near_wall_dU_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*Obj.ks);
                Near_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro>=0 & WT10VorX.z_pro <= m*Obj.ks);

                Far_wall_dU_pro_WT10 = WT10VorX.dU_pro (WT10VorX.z_pro > m*Obj.ks);
                Far_wall_d_pro_WT10 = WT10VorX.d_pro (WT10VorX.z_pro > m*Obj.ks);


                R1 = corrcoef(0.5*Near_wall_d_pro_WT10,-Near_wall_dU_pro_WT10);
                rho_near = abs(R1(1,2)); % Target (Spearman) correlation %

                R2 = corrcoef(0.5*Far_wall_d_pro_WT10,-Far_wall_dU_pro_WT10);
                rho_Far = abs(R2(1,2)); % Target (Spearman) correlation %
            end
            delUwOu_tau_near = load('delUwOu_tau_near.mat').delUwOu_tau_near;
            delUwOu_tau_Far = load('delUwOu_tau_Far.mat').delUwOu_tau_Far;
            rwOlambda_T_near = load('rwOlambda_T_near.mat').rwOlambda_T_near;
            rwOlambda_T_Far = load('rwOlambda_T_Far.mat').rwOlambda_T_Far;
            ii = 1;
            iii = 1;

            z_prog = zeros(1,1);
            z_retro = zeros(1,1);
            z_prog_near = zeros(1,1);
            pr = 1;
            rt = 1;
            pr_near = 1;
            rho_Near_wall = load('samples_Near.mat').samples_Near;
            rho_Far_wall = load('samples_Far.mat').samples_Far;
            xi = 1;
            xii = 1;
            Dist_close_wall_VorX = zeros(1);
            ind_close_wall_VorX = 1;
            for S=1:size(Obj.HRVFu,3)

                xend_near_wall_prev = 0;
                temp_ind_close_wall_VorX = 1;
                while xend_near_wall_prev < max(Obj.HRVFx)% Attached wall vortices

                    D1 = 1;
                    D2 = 1;
                    N = 1;
                    fac = 1.0;% should be more than 1.0
                    r_omega = rwOlambda_T_near(ii)*Obj.Delx; %choose from near wall distribution
                    zc_near_wall = N*r_omega + Obj.HRVFz(1);
                    xc_near_wall = N*r_omega + xend_near_wall_prev;
                    xend_near_wall_prev = xc_near_wall + fac*N*r_omega;
                    [indzmid] = find(Obj.HRVFz>=zc_near_wall,1);
                    [indxmid] = find(Obj.HRVFx>=xc_near_wall,1);


                    %                     logic = randi([0, 1]);

                    %                     if Obj.HRVFuprime(indzmid,indxmid,S)<0
                    %                        ii = ii + 1;
                    %                        xi= xi+1;
                    %                        xii = xii+1;
                    %                        continue
                    %                     end

                    if randi([0, 1])==0
                        ii = ii + 1;
                        xi= xi+1;
                        xii = xii+1;
                        continue
                    end

                    if temp_ind_close_wall_VorX > 1
                        Dist_close_wall_VorX(ind_close_wall_VorX-1) = xc_near_wall-prev_pos_xc_near_wall_VorX;
                    end
                    prev_pos_xc_near_wall_VorX = xc_near_wall;
                    temp_ind_close_wall_VorX = temp_ind_close_wall_VorX +1;
                    ind_close_wall_VorX = ind_close_wall_VorX +1;

                    if zc_near_wall<= 2*Obj.ks

                        Trans = [1, rho_Near_wall(xi); 0, sqrt(1 - rho_Near_wall(xi)^2)];
                        xi= xi+1;
                    else

                        Trans = [1, rho_Far_wall(xii); 0, sqrt(1 - rho_Far_wall(xii)^2)];
                        xii = xii+1;
                    end

                    Gama=2*pi*delUwOu_tau_near(ii)*Obj.u_tau*r_omega/(1-exp(-1));  %choose from near wall distribution
                    ii = ii + 1;
                    [indz] = find(Obj.HRVFz>=zc_near_wall-N*r_omega &...
                        Obj.HRVFz<=zc_near_wall+N*r_omega);
                    [indx] = find(Obj.HRVFx>=xc_near_wall-N*r_omega &...
                        Obj.HRVFx<=xc_near_wall+N*r_omega);
                    [X,Z] = meshgrid(Obj.HRVFx(indx),Obj.HRVFz(indz));
                    r = sqrt((X-xc_near_wall).^2+(Z-zc_near_wall).^2);
                    %                     r(r>N*1.01*r_omega)=0;
                    z_prog_near(pr_near,1) = zc_near_wall;
                    pr_near = pr_near +1;

                    uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
                    nan_logical_array = isnan(uazi);
                    [nan_indices] = find(nan_logical_array);
                    uazi(nan_indices) = 0;

                    Theta = atan2((Z-zc_near_wall),(X-xc_near_wall));
                    nan_logical_array = isnan(Theta);
                    [nan_indices] = find(nan_logical_array);
                    Theta(nan_indices) = 0;

                    mask1 = (Theta >= 0 & Theta < pi/2 & r~=0) ;

                    mask2 =(Theta >= pi/2 & Theta <= pi+0.1 & r~=0);

                    mask3 =((Theta >= -pi & Theta <= -pi/2 & r~=0));

                    mask4 = (Theta >= -pi/2 & Theta <= 0 & r~=0) ;


                    u_vor = uazi.*sin(Theta);
                    u_vor(uazi~=0)=u_vor(uazi~=0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                    w_vor = -uazi.*cos(Theta);
                    %                     w_vor(uazi>0)=w_vor(uazi>0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                    A1 = zeros(size(w_vor(mask1),1),2);
                    A2 = zeros(size(w_vor(mask2),1),2);
                    A3 = zeros(size(w_vor(mask3),1),2);
                    A4 = zeros(size(w_vor(mask4),1),2);
                    A1(:,1) = u_vor(mask1);
                    A1(:,2) = w_vor(mask1);
                    A2(:,1) = u_vor(mask2);
                    A2(:,2) = w_vor(mask2);
                    A3(:,1) = u_vor(mask3);
                    A3(:,2) = w_vor(mask3);
                    A4(:,1) = u_vor(mask4);
                    A4(:,2) = w_vor(mask4);


                    B1 = A1* Trans;
                    B2 = A2* Trans;
                    B3 = A3* Trans;
                    B4 = A4* Trans;
                    u_vor_trans = zeros(size(u_vor));
                    w_vor_trans = zeros(size(w_vor));
                    u_vor_trans(mask1) = B1(:,1);
                    w_vor_trans(mask1) = B1(:,2);

                    u_vor_trans(mask2) = B2(:,1);
                    w_vor_trans(mask2) = B2(:,2);

                    u_vor_trans(mask3) = B3(:,1);
                    w_vor_trans(mask3) = B3(:,2);

                    u_vor_trans(mask4) = B4(:,1);
                    w_vor_trans(mask4) = B4(:,2);


                    if isempty(u_vor)
                        continue
                    end

                    Obj.HRVFu(indz,indx,S) = Obj.HRVFu(indz,indx,S)+ D1*u_vor_trans;
                    %                     Obj.HRVFuprime(indz,indx,S) = Obj.HRVFuprime(indz,indx,S)+ u_vor_trans;
                    Obj.HRVFw(indz,indx,S) = Obj.HRVFw(indz,indx,S)+ D2*w_vor_trans;

                    %
                    %                     Obj.HRVFu(indz,indx,S) = Obj.HRVFu(indz,indx,S)+ D1*u_vor;
                    % %                     Obj.HRVFuprime(indz,indx,S) = Obj.HRVFuprime(indz,indx,S)+ u_vor;
                    %                     Obj.HRVFw(indz,indx,S) = Obj.HRVFw(indz,indx,S)+ D2*w_vor;
                end
                for p=1:size(Obj.x_shear,1)
                    if isempty(Obj.x_shear{p,S})
                        continue
                    end
                    for i =1:length(Obj.x_shear{p,S})
                        xcs = Obj.x_shear{p,S}{1,i};

                        zcs = Obj.z_shear{p,S}{1,i};

                        rot = -sign(Obj.Delta_u_m{p,S}{1,i});
                        N1 = 1;%2This should be 2 otherwise E1, E2 has no significant effect
                        E1 = 1;
                        E2 = 1;
                        [Obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                            lambdacivorX(Obj,1,ii,iii,S,xcs,zcs,rot,rwOlambda_T_near,rwOlambda_T_Far,...
                            abs(Obj.Delta_u_m{p,S}{1,i})*0.5/Obj.u_tau,abs(Obj.Delta_u_m{p,S}{1,i})*0.5/Obj.u_tau,...
                            z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                            rho_Far_wall(xii), xii, rho_near, rho_Far);
                    end

                end
            end
            Obj.HRVFuprime = Obj.HRVFu - mean(Obj.HRVFu,3);
            %             Obj.HRVFuprime = Obj.HRVFu - Obj.u_tau/kappa*log(Obj.HRVFz/Obj.z0);
            numb_prog = zeros(size(Obj.HRVFz)-1);
            numb_prog_near = zeros(size(Obj.HRVFz)-1);
            numb_retro = zeros(size(Obj.HRVFz)-1);
            for k = 1:size(Obj.HRVFz,1)-1
                numb_prog(k,1) = sum(z_prog>=Obj.HRVFz(k,1)&...
                    z_prog<=Obj.HRVFz(k+1,1),1);
                numb_prog_near(k,1) = sum(z_prog_near>=Obj.HRVFz(k,1)&...
                    z_prog_near<=Obj.HRVFz(k+1,1),1);
                numb_retro(k,1) = sum(z_retro>=Obj.HRVFz(k,1)&...
                    z_retro<=Obj.HRVFz(k+1,1),1);
            end
            figure
            subplot(1,3,1)
            plot(numb_prog/(Obj.HRVFx(end)*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.delta))
            set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
                'YGrid','on')
            ylabel('z/$\delta$','Interpreter','Latex','FontSize',14);
            xlabel('Number of prograde/(Area)','Interpreter','Latex','FontSize',14);
            subplot(1,3,2)
            plot(numb_prog_near/(Obj.HRVFx(end)*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.delta))
            set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
                'YGrid','on')
            ylabel('z/$\delta$','Interpreter','Latex','FontSize',14);
            xlabel('Number of prograde/(Area)','Interpreter','Latex','FontSize',14);
            subplot(1,3,3)
            plot(numb_retro/(Obj.HRVFx(end)*size(Obj.HRVFu,3)),...
                (Obj.HRVFz(1:end-1) + Obj.HRVFz(2:end)) / (2*Obj.delta),'r')
            set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
                'YGrid','on')
            ylabel('z/$\delta$','Interpreter','Latex','FontSize',14);
            xlabel('Number of retrograde/(Area)','Interpreter','Latex','FontSize',14);
            disp(['The mean value of vortex distance normalized with ks is: ', num2str(mean(Dist_close_wall_VorX,2)/(Obj.ks))]);
            disp(['The std value of vortex distance normalized with ks is: ', num2str(std(Dist_close_wall_VorX,0,2)/(Obj.ks))]);
        end
        function Spectral_analysis(Obj)%%%DO SPECTRA FOR H-W
            if Obj.type == "Exp"
                if Obj.name =="SonicASL"
                    win12del=floor(12*Obj.delta./mean(Obj.u,1)*Obj.fs);
                    Obj.wavenumb = cell(size(win12del,1),1);%[];%zeros(size(win12del,1));
                    Obj.PowSpecDen_k = cell(size(win12del,1),1);%[];%zeros(size(win12del,1),1);
                    Obj.epsilon_spec = zeros(size(win12del,1),1);
                    Obj.eta_spec = zeros(size(win12del,1),1);
                    for i =1:size(win12del,1)
                        window=hann(win12del(i,1));
                        NFFT=2^(nextpow2(win12del(i,1)));


                        [PSD_f,freq]=pwelch(Obj.uprime(:,1),window,NFFT/2,win12del(i,1),Obj.fs); %PSD_f[m^2/s]

                        k = 2*pi*freq/mean(Obj.u(:,1),1);%
                        Obj.wavenumb{i,1} = k;

                        Akfr = trapz(freq,PSD_f);%%[m^2/s^3]=[m^2/s^2]*[1/s]
                        Q =  mean(Obj.uprime(:,1).*Obj.uprime(:,1),1)/Akfr;%%[s]=[m^2/s^2]/[m^2/s^3]
                        PSD_f = PSD_f * Q;%%[m^2/s]=[m^2/s^2]*[s]
                        Obj.PowSpecDen_k{i,1} = PSD_f*mean(Obj.u(:,1),1)/(2*pi);%%[m^3/s^2]=[m^2/s]*[m/s]
                        Obj.epsilon_spec(i,1) = 15*Obj.nu*trapz(k,(k.^2).*Obj.PowSpecDen_k{i,1});
                        Obj.eta_spec(i,1) = (Obj.nu^3/Obj.epsilon_spec(i,1))^0.25;
                    end
                else
                    win12del=floor(12*Obj.delta./mean(Obj.u,2)*Obj.fs);
                    Obj.wavenumb = cell(size(win12del,1),1);%[];%zeros(size(win12del,1));
                    Obj.PowSpecDen_k = cell(size(win12del,1),1);%[];%zeros(size(win12del,1),1);
                    Obj.epsilon_spec = zeros(size(win12del,1),1);
                    Obj.eta_spec = zeros(size(win12del,1),1);
                    for i =1:size(win12del,1)
                        window=hann(win12del(i,1));
                        NFFT=2^(nextpow2(win12del(i,1)));


                        [PSD_f,freq]=pwelch(Obj.uprime(i,:),window,NFFT/2,win12del(i,1),Obj.fs); %PSD_f[m^2/s]

                        k = 2*pi*freq/mean(Obj.u(i,:),2);%
                        Obj.wavenumb{i,1} = k;

                        Akfr = trapz(freq,PSD_f);%%[m^2/s^3]=[m^2/s^2]*[1/s]
                        Q =  mean(Obj.uprime(i,:).*Obj.uprime(i,:),2)/Akfr;%%[s]=[m^2/s^2]/[m^2/s^3]
                        PSD_f = PSD_f * Q;%%[m^2/s]=[m^2/s^2]*[s]
                        Obj.PowSpecDen_k{i,1} = PSD_f*mean(Obj.u(i,:),2)/(2*pi);%%[m^3/s^2]=[m^2/s]*[m/s]
                        Obj.epsilon_spec(i,1) = 15*Obj.nu*trapz(k,(k.^2).*Obj.PowSpecDen_k{i,1});
                        Obj.eta_spec(i,1) = (Obj.nu^3/Obj.epsilon_spec(i,1))^0.25;
                    end
                end
            else
                concat_uprime = reshape(Obj.HRVFuprime,size(Obj.HRVFuprime, 1), []);
                ksmax=floor(2*pi/(Obj.HRVFx(2)-Obj.HRVFx(1)));

                win12del=floor(12*Obj.delta/(Obj.HRVFx(2)-Obj.HRVFx(1)));
                window=hann(win12del);
                NFFT=2^(nextpow2(win12del));
                [PSD_k,k]=pwelch(concat_uprime',window,NFFT/2,win12del,ksmax); %PSD_f[m^2/s]
                Apk = trapz(k,PSD_k);
                Q = mean(var(Obj.HRVFu,0,3),2)./Apk';  %%[m^3/s^2]=[m^2/s]*[m^2/s^2]*1/[m/s]
                Obj.PowSpecDen_k=(PSD_k.*Q')';
                %                 trapz(k,PSD_k(:,10)*Q(10))

                % PSD_kWT7 = 2*PSD_kWT7(1:ceil(w/2));
                % kWT7 = kWT7(1:ceil(w/2));
                Obj.wavenumb = k;
                Obj.epsilon_spec = 15*Obj.nu*trapz(k,PSD_k.*(k.^2));
                Obj.eta_spec = (Obj.nu^3./Obj.epsilon_spec).^0.25;

                %                 produWT7=mean(uprimefieldGenWT7(ptsWTGen(1,4),:,fN).*wfieldGenWT7(ptsWTGen(1,4),:,fN),2)*...
                %                     (mean(ufieldGenWT7(ptsWTGen(1,4),:,fN),2)-mean(ufieldGenWT7(ptsWTGen(1,4)-1,:,fN),2))/...
                %                     (meanzWT7(ptsWTGen(1,4),1)-meanzWT7(ptsWTGen(1,4)-1,1));
            end

        end
        function Structure_function(Obj)
            if Obj.type == "Exp"
                if Obj.name == "SonicASL"
                    Obj.D11=cell(size(Obj.u,2),1);
                    Obj.r11=cell(size(Obj.u,2),1);
                    Obj.epsilon_str = cell(size(Obj.u,2),1);
                    Obj.eta_str = cell(size(Obj.u,2),1);
                    for elev = 1:size(Obj.z,1)
                        delx = mean(Obj.u(:,1),1)*1/Obj.fs;
                        win12del = ceil(12*Obj.delta/(delx));
                        Numb_frame = floor(size(Obj.u,1)/win12del);
                        Obj.D11{elev,1}=cell(1,Numb_frame);
                        Obj.r11{elev,1}=cell(1,1);
                        Obj.epsilon_str{elev,1} = cell(1,Numb_frame);
                        Obj.eta_str{elev,1} = cell(1,Numb_frame);
                        for r=0:win12del
                            Obj.r11{elev,1}{r+1}=(delx)*r;
                        end
                        Obj.r11{elev,1} = [Obj.r11{elev,1}{:}];
                        for S = 1:Numb_frame
                            for r = 0:win12del
                                Obj.D11{elev,1}{1,S}{r+1}=mean((Obj.u((S-1)*win12del+1:S*win12del+1-r,1)-...
                                    Obj.u((S-1)*win12del+1+r:S*win12del+1,1)).^2,1);
                            end
                            Obj.D11{elev,1}{1,S} = [Obj.D11{elev,1}{1,S}{:}];

                            [max_val]=max(Obj.D11{elev,1}{1,S}.*Obj.r11{elev,1}.^(-2/3),[],2);
                            Obj.epsilon_str{elev,1}{1,S} = (2./max_val).^(-3/2);
                            Obj.eta_str{elev,1}{1,S} = (Obj.nu^3./Obj.epsilon_str{elev,1}{1,S}).^0.25;
                        end
                    end
                else
                    Obj.D11=cell(size(Obj.u,1),1);
                    Obj.r11=cell(size(Obj.u,1),1);
                    Obj.epsilon_str = cell(size(Obj.u,1),1);
                    Obj.eta_str = cell(size(Obj.u,1),1);
                    for elev = 1:size(Obj.z,1)
                        delx = mean(Obj.u(elev,:),2)*1/Obj.fs;
                        win12del = ceil(12*Obj.delta/(delx));
                        Numb_frame = floor(size(Obj.u,2)/win12del);
                        Obj.D11{elev,1}=cell(1,Numb_frame);
                        Obj.r11{elev,1}=cell(1,1);
                        Obj.epsilon_str{elev,1} = cell(1,Numb_frame);
                        Obj.eta_str{elev,1} = cell(1,Numb_frame);
                        for r=0:win12del
                            Obj.r11{elev,1}{r+1}=(delx)*r;
                        end
                        Obj.r11{elev,1} = [Obj.r11{elev,1}{:}];
                        for S = 1:Numb_frame
                            for r = 0:win12del
                                Obj.D11{elev,1}{1,S}{r+1}=mean((Obj.u(elev,(S-1)*win12del+1:S*win12del+1-r)-...
                                    Obj.u(elev,(S-1)*win12del+1+r:S*win12del+1)).^2,2);
                            end
                            Obj.D11{elev,1}{1,S} = [Obj.D11{elev,1}{1,S}{:}];

                            [max_val]=max(Obj.D11{elev,1}{1,S}.*Obj.r11{elev,1}.^(-2/3),[],2);
                            Obj.epsilon_str{elev,1}{1,S} = (2./max_val).^(-3/2);
                            Obj.eta_str{elev,1}{1,S} = (Obj.nu^3./Obj.epsilon_str{elev,1}{1,S}).^0.25;
                        end
                    end
                end
            else
                concat_uprime = reshape(Obj.HRVFuprime,size(Obj.HRVFuprime, 1), []);
                win12del = ceil(12*Obj.delta/(Obj.HRVFx(2)-Obj.HRVFx(1)));
                Numb_frame = floor(size(concat_uprime,2)/win12del);
                Obj.D11=zeros(size(concat_uprime,1),win12del+1,Numb_frame);
                Obj.r11=zeros(1,win12del+1);
                Obj.epsilon_str = zeros(size(concat_uprime,1),Numb_frame);
                Obj.eta_str = zeros(size(concat_uprime,1),Numb_frame);
                for r=0:win12del
                    Obj.r11(1,r+1)=(Obj.HRVFx(2)-Obj.HRVFx(1))*r;
                end
                if Obj.name == "GenWT7" % Since it is time consuming, I did for the elevations that I want to have plot
                    Z = [4,203];
                elseif Obj.name == "GenWT10"
                    Z = [21,200];
                else
                    Z = [6,64];
                end
                for S = 1:Numb_frame
                    for elev = Z%1:size(Obj.z,1)
                        for r=0:size(Obj.r11,2)-1
                            Obj.D11(elev,r+1,S)=mean((concat_uprime(elev,(S-1)*win12del+1:S*win12del+1-r)-...
                                concat_uprime(elev,(S-1)*win12del+1+r:S*win12del+1)).^2,2);
                        end
                    end
                    [max_val]=max(Obj.D11(:,:,S).*Obj.r11.^(-2/3),[],2);
                    Obj.epsilon_str(:,S) = (2./max_val).^(-3/2);
                    Obj.eta_str(:,S)=(Obj.nu^3./Obj.epsilon_str(:,S)).^0.25;
                end
            end

        end
        function Lambda_T_computation(Obj, Aux_Obj)
            if Obj.type == "Exp"
                if Obj.name ~= "SonicASL"
                    Obj.Lambda_T = zeros(size(Obj.z, 1),1);
    
                    for i = 1:size(Obj.z, 1)                    
                        Obj.Lambda_T(i) = rms(Obj.uprime(i,:),2)*sqrt(15*Obj.nu/mean([Obj.epsilon_str{i, :}{:}],2));
                    end
                else
                    Obj.Lambda_T = rms(Obj.uprime(:,1),1)*sqrt(15*Obj.nu/mean([Obj.epsilon_str{1, :}{:}],2));
                end
            else
                if Obj.name ~= "GenASL"               
                    Obj.Lambda_T = interp1(Aux_Obj.z, Aux_Obj.Lambda_T, Obj.z, 'linear', 'extrap');
                else
                    Obj.Lambda_T = Aux_Obj.Lambda_T;
                end

            end

        end
        function Newfunction(Obj, propName)
            if ~isprop(Obj, propName)
                addprop(Obj, propName);
            end
            Obj.(propName) = zeros(size(Obj.u));

        end

    end
end