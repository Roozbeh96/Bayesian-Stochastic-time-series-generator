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
%                 if isnan(factor)
%                     Obj.factor = floor(DelxGen/DelzGen);
%                 else
%                     Obj.factor = factor;
%                 end
                %                 Obj.HRVFu=zeros(size(Obj.u,1), Obj.factor*(L) +1, Numrator);
                %                 Obj.HRVFw=zeros(size(Obj.u,1), Obj.factor*(L) +1, Numrator);
                %                 Obj.HRVFx= 0:DelxGen/Obj.factor: (L)*DelxGen;
                %                 Obj.HRVFz=Obj.z;
                xorg = 0:DelxGen:(L)*DelxGen;
                Obj.x = xorg;
                %                 xi = Obj.HRVFx;
                xi = size(Obj.HRVFu,2)*size(Obj.HRVFu,3);
                res = (Obj.HRVFx(2)-Obj.HRVFx(1));

%                 size_minmaxvec = 1.0e6;
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
%                 save('minmaxvec_Long_Long.mat','minmaxvec')
%                 save('minmaxvec_Long.mat','minmaxvec')
%                 data = load('minmaxvec_short.mat');
%                 data = load('minmaxvec.mat');
                data = load('minmaxvec_Long.mat');
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

    end
end