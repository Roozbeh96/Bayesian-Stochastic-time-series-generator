function [Obj,ii,z_prog, z_retro, pr, rt] = lambdacivorX_mac(Obj,i,ii,S,xcs,zcs,rot,dwOlambda_T,...
                        delUwOu_tau, z_prog, z_retro, pr, rt, Trans, N, C1, C2)


r_omega = dwOlambda_T(ii)*Obj.Delx*0.5;
if Obj.HRVFz(1)>=zcs(i)-N*r_omega
    r_omega = zcs(i)-Obj.HRVFz(1);
end

if Obj.HRVFz(end)<=zcs(i)+N*r_omega
    r_omega = Obj.HRVFz(end)-zcs(i);
end

% value of the N should be calculated based on the
% r_omega(it should be dissipated) but here is random
% between 1 to 10
Gama=2*pi*0.5*delUwOu_tau(ii)*Obj.u_tau*r_omega;
ii = ii + 1;
[indz] = find(Obj.HRVFz>=zcs(i)-N*r_omega &...
    Obj.HRVFz<=zcs(i)+N*r_omega);
[indx] = find(Obj.HRVFx>=xcs(i)-N*r_omega &...
    Obj.HRVFx<=xcs(i)+N*r_omega);
[X,Z] = meshgrid(Obj.HRVFx(indx),Obj.HRVFz(indz));
r = sqrt((X-xcs(i)).^2+(Z-zcs(i)).^2);
r(r>N*r_omega)=0;



uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
nan_logical_array = isnan(uazi);
[nan_indices] = find(nan_logical_array);
uazi(nan_indices) = 0;

Theta = atan2((Z-zcs(i)),(X-xcs(i)));
nan_logical_array = isnan(Theta);
[nan_indices] = find(nan_logical_array);
Theta(nan_indices) = 0;

mask1 = (Theta >= 0 & Theta < pi/2 & r~=0) ;

mask2 =(Theta >= pi/2 & Theta <= pi+0.1 & r~=0);

mask3 =((Theta >= -pi & Theta <= -pi/2 & r~=0));

mask4 = (Theta >= -pi/2 & Theta <= 0 & r~=0) ;

%                     r_normalized = r / max(r(:));
%                     cmap1 = getPyPlot_cMap('Reds_r', 8,[],...
%                           '"C:\Users\ehsan010\Anaconda3\envs\General\python.exe"');
%                     cmap1(end-1:end,:)=[];
%                     cmap2 = getPyPlot_cMap('Greens_r', 8, [],...
%                           '"C:\Users\ehsan010\Anaconda3\envs\General\python.exe"');
%                     cmap2(end-1:end,:)=[];
%                     cmap3 = getPyPlot_cMap('Blues_r', 8, [],...
%                           '"C:\Users\ehsan010\Anaconda3\envs\General\python.exe"');
%                     cmap3(end-1:end,:)=[];
%                     cmap4 = getPyPlot_cMap('Wistia_r', 6, [],...
%                         '"C:\Users\ehsan010\Anaconda3\envs\General\python.exe"');
%                     color_indices1 = round((1-r_normalized) * (size(cmap1, 1)/1 - 1)) + 1;
%                     color_indices2 = round((1-r_normalized) * (size(cmap2, 1)/1 - 1)) + 1;
%                     color_indices3 = round((1-r_normalized) * (size(cmap3, 1)/1 - 1)) + 1;
%                     color_indices4 = round((1-r_normalized) * (size(cmap4, 1)/1 - 1)) + 1;
%
%                     figure
%                     set(gcf,'Position',[1291,520,540,450])
%
%                     scatter((X(mask1)-xcs(i))/r_omega, (Z(mask1)-zcs(i))/r_omega...
%                         , 10, cmap1(color_indices1(mask1), :), 'filled');
%                     hold on
%                     scatter((X(mask2)-xcs(i))/r_omega, (Z(mask2)-zcs(i))/r_omega,...
%                         10, cmap2(color_indices2(mask2), :), 'filled');
%                     scatter((X(mask3)-xcs(i))/r_omega, (Z(mask3)-zcs(i))/r_omega,...
%                         10, cmap3(color_indices3(mask3), :), 'filled');
%                     scatter((X(mask4)-xcs(i))/r_omega, (Z(mask4)-zcs(i))/r_omega,...
%                         10, cmap4(color_indices4(mask4), :), 'filled');
%                     set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
%                         'YGrid','on')
%                     xlabel('x/$r_{\omega}$','Interpreter','Latex','FontSize',14);
%                     ylabel('z/$r_{\omega}$','Interpreter','Latex','FontSize',14);
%                     axis equal
%                     xlim([-2.2 2.2])
%                     ylim([-2.2 2.2])
%                     box on
%                     Theta_t = pi / 10;
%                     mask1 = (Theta > -Theta_t & Theta < Theta_t) ;
%
%                     mask2 =(Theta > pi-Theta_t);
%
%                     mask3 =(Theta > pi/2-Theta_t & Theta < pi/2+Theta_t);
%
%                     mask4 = (Theta > -pi/2-Theta_t & Theta < -pi/2+Theta_t) ;
%
%                     mask5 = (Theta < -pi+Theta_t);
%
%                     % being prograde or retrograde also can be stochastic
%                     % for prograde
% %                     prograd = logical(binornd(1, 0.5));
%                     uazi_aug = zeros(size(Theta));
%                     uazi_aug(mask1) = abs(2 * cos(pi * Theta(mask1) / (2 * Theta_t)));
%                     uazi_aug(mask2) = abs(2 * cos(pi * (Theta(mask2)-pi) / (2 * Theta_t)));
%                     uazi_aug(mask3) = abs(2 * cos(pi * (Theta(mask3)-pi/2) / (2 * Theta_t)));
%                     uazi_aug(mask4) = abs(2 * cos(pi * (Theta(mask4)+pi/2) / (2 * Theta_t)));
%                     uazi_aug(mask5) = abs(2 * cos(pi * (Theta(mask5)+pi) / (2 * Theta_t)));


if rot(i)<0
    z_prog(pr,1) = zcs(i);
    pr = pr +1;
    u_vor = uazi.*sin(Theta);
    w_vor = -uazi.*cos(Theta);
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
    
    %                         figure
    %                         set(gcf,'Position',[1291,520,540,450])
    %                         scatter(w_vor(mask1)/Obj.u_tau, u_vor(mask1)/Obj.u_tau...
    %                             , 10, cmap1(color_indices1(mask1), :), 'filled');
    %                         hold on
    %                         scatter(w_vor(mask2)/Obj.u_tau, u_vor(mask2)/Obj.u_tau...
    %                             , 10, cmap2(color_indices2(mask2), :), 'filled');
    %                         scatter(w_vor(mask3)/Obj.u_tau, u_vor(mask3)/Obj.u_tau...
    %                             , 10, cmap3(color_indices3(mask3), :), 'filled');
    %                         scatter(w_vor(mask4)/Obj.u_tau, u_vor(mask4)/Obj.u_tau...
    %                             , 10, cmap4(color_indices4(mask4), :), 'filled');
    %                         set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
    %                             'YGrid','on')
    %                         xlabel('w/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         ylabel('u/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         axis equal
    %                         xlim([-0.8 0.8])
    %                         ylim([-0.8 0.8])
    %                         box on
    B1 = A1* Trans;
    B2 = A2* Trans;
    B3 = A3* Trans;
    B4 = A4* Trans;
    u_vor_trans = zeros(size(u_vor));
    w_vor_trans = zeros(size(w_vor));
    u_vor_trans(mask1) = B1(:,1);%%any coefficient should be multiplied here
    w_vor_trans(mask1) = B1(:,2);
    
    u_vor_trans(mask2) = B2(:,1);
    w_vor_trans(mask2) = B2(:,2);
    
    u_vor_trans(mask3) = B3(:,1);
    w_vor_trans(mask3) = B3(:,2);
    
    u_vor_trans(mask4) = B4(:,1);
    w_vor_trans(mask4) = B4(:,2);
    
    %                         figure
    %                         set(gcf,'Position',[1291,520,540,450])
    %                         scatter(B1(:,2)/Obj.u_tau, B1(:,1)/Obj.u_tau...
    %                             , 10, cmap1(color_indices1(mask1), :), 'filled');
    %                         hold on
    %                         scatter(B2(:,2)/Obj.u_tau, B2(:,1)/Obj.u_tau...
    %                             , 10, cmap2(color_indices2(mask2), :), 'filled');
    %                         scatter(B3(:,2)/Obj.u_tau, B3(:,1)/Obj.u_tau...
    %                             , 10, cmap3(color_indices3(mask3), :), 'filled');
    %                         scatter(B4(:,2)/Obj.u_tau, B4(:,1)/Obj.u_tau...
    %                             , 10, cmap4(color_indices4(mask4), :), 'filled');
    %                         set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
    %                             'YGrid','on')
    %                         xlabel('w/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         ylabel('u/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         axis equal
    %                         xlim([-0.8 0.8])
    %                         ylim([-0.8 0.8])
    %                         box on
    
    %                         figure
    %                         Gammaoseen = sqrt(u_vor.^2+w_vor.^2)*2*pi.*r./(1-exp(-r.^2./r_omega^2));
    %                         Gammatrans = sqrt(u_vor_trans.^2+w_vor_trans.^2)*2*pi.*r./(1-exp(-r.^2./r_omega^2));
    %                         contourf((X-xcs(i))./r_omega,...
    %                             (Z-zcs(i))./r_omega,...
    %                             Gammatrans-Gammaoseen...
    %                             ,50,'LineStyle','none',...
    %                             'LineColor','none','EdgeColor','none');
    %                         % cmap2 = plasma(100);
    % %                         cmap2 = crameri('turku',100);
    %                         cmap2 = getPyPlot_cMap('RdGy_r', [], [],...
    %                             '"C:\Users\ehsan010\Anaconda3\envs\General\python.exe"');
    % %                         cmap2 = flip(cmap2);
    %                         colormap(gca,cmap2)
    %                         % C2=cmap2;
    %                         % C3=[0,0,0;C2;1,1,1];
    %                         caxis([-0.003 0.003])
    %                         % C3 = flip(C3);
    %                         % colormap(gca,C3)
    %                         hcb2=colorbar;
    %                         % hcb2.Position=[0.954656862745097,0.407407407407417,0.01593137254902,0.229885057471265];
    %                         title(hcb2,'$\Gamma(Trans)$-$\Gamma(Oseen)$','Interpreter','Latex','FontSize',15)...,'Position',...
    %                         %     [4.999994333833513,142.5300008152491,0]);
    %                         hcb2.TickLabelInterpreter = 'latex';
    %                         hcb2.TickLength = 0;
    %                         hcb2.Label.Interpreter = 'latex';
    %                         hcb2.AxisLocation='out';
    %                         set(gca,'TickLabelInterpreter','latex','FontSize',13,...
    %                             'XGrid','on','YGrid','on')
    %                         xlabel('x/r$_{\omega}$','Interpreter','Latex','FontSize',14);
    %                         ylabel('z/r$_{\omega}$','Interpreter','Latex','FontSize',14);
    %                         axis equal
    %                         %xlim([-2.2 2.2])
    %                         %ylim([-2.2 2.2])
    %
    %                                                 figure
    %                                                 quiver((X(1:2:end,1:2:end)-xcs(i))./r_omega,...
    %                                                     (Z(1:2:end,1:2:end)-zcs(i))./r_omega,...
    %                                                     u_vor_trans(1:2:end,1:2:end),w_vor_trans(1:2:end,1:2:end)...
    %                                                     ,2,'color',[1.00,0.00,0.00],...
    %                                                     'MaxHeadSize',0.2,'LineWidth',1.5);
    %                                                 set(gca,'TickLabelInterpreter','latex','FontSize',13,...
    %                                                     'XGrid','on','YGrid','on')
    %                                                 xlabel('x/r$_{\omega}$','Interpreter','Latex','FontSize',14);
    %                                                 ylabel('z/r$_{\omega}$','Interpreter','Latex','FontSize',14);
    %                                                 axis equal
    %                                                 xlim([-2.2 2.2])
    %                                                 ylim([-2.2 2.2])
    
else
    z_retro(rt,1) = zcs(i);
    rt = rt+1;
    u_vor = -uazi.*sin(Theta);
    w_vor = uazi.*cos(Theta);
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
    
    %                         figure
    %                         set(gcf,'Position',[1291,520,540,450])
    %                         scatter(w_vor(mask1)/Obj.u_tau, u_vor(mask1)/Obj.u_tau...
    %                             , 10, cmap1(color_indices1(mask1), :), 'filled');
    %                         hold on
    %                         scatter(w_vor(mask2)/Obj.u_tau, u_vor(mask2)/Obj.u_tau...
    %                             , 10, cmap2(color_indices2(mask2), :), 'filled');
    %                         scatter(w_vor(mask3)/Obj.u_tau, u_vor(mask3)/Obj.u_tau...
    %                             , 10, cmap3(color_indices3(mask3), :), 'filled');
    %                         scatter(w_vor(mask4)/Obj.u_tau, u_vor(mask4)/Obj.u_tau...
    %                             , 10, cmap4(color_indices4(mask4), :), 'filled');
    %                         set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
    %                             'YGrid','on')
    %                         xlabel('w/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         ylabel('u/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         axis equal
    %                         xlim([-0.8 0.8])
    %                         ylim([-0.8 0.8])
    %                         box on
    
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
    
    %                         figure
    %                         set(gcf,'Position',[1291,520,540,450])
    %                         scatter(B1(:,2)/Obj.u_tau, B1(:,1)/Obj.u_tau...
    %                             , 10, cmap1(color_indices1(mask1), :), 'filled');
    %                         hold on
    %                         scatter(B2(:,2)/Obj.u_tau, B2(:,1)/Obj.u_tau...
    %                             , 10, cmap2(color_indices2(mask2), :), 'filled');
    %                         scatter(B3(:,2)/Obj.u_tau, B3(:,1)/Obj.u_tau...
    %                             , 10, cmap3(color_indices3(mask3), :), 'filled');
    %                         scatter(B4(:,2)/Obj.u_tau, B4(:,1)/Obj.u_tau...
    %                             , 10, cmap4(color_indices4(mask4), :), 'filled');
    %                         set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on',...
    %                             'YGrid','on')
    %                         xlabel('w/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         ylabel('u/$u_{\tau}$','Interpreter','Latex','FontSize',14);
    %                         axis equal
    %                         xlim([-0.8 0.8])
    %                         ylim([-0.8 0.8])
    %                         box on
end
Obj.HRVFu(indz,indx,S) = Obj.HRVFu(indz,indx,S)+ C1*u_vor_trans;
%                     Obj.HRVFuprime(indz,indx,S) = Obj.HRVFuprime(indz,indx,S)+ u_vor_trans;
Obj.HRVFw(indz,indx,S) = Obj.HRVFw(indz,indx,S)+ C2*w_vor_trans;

end

