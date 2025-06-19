
N = 4;
dwOlambda_T = 0.4;
lambda_T = 0.01;
delUwOu_tau = 2;
u_tau = 0.4;
r_omega = dwOlambda_T*lambda_T*0.5;
Gama=2*pi*0.5*delUwOu_tau*u_tau*r_omega;


[X,Z] = meshgrid(-N*r_omega:0.0001:N*r_omega,-N*r_omega:0.0001:N*r_omega);
r = sqrt((X).^2+(Z).^2);
U = 0.00;
m = U/(2*N*r_omega);
U_field = m*(Z+N*r_omega);
% U_field = U*ones(size(X));

% U_field(X>r_omega|X<-r_omega)=0;
% U_field(Z>r_omega|Z<-r_omega)=0;
U_field(r>N*r_omega)=0;
r(r>N*r_omega)=0;

ome = 1;
uazi = r.*ome;
% uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
nan_logical_array = isnan(uazi);
[nan_indices] = find(nan_logical_array);
uazi(nan_indices) = 0;

Theta = atan2((Z),(X));
nan_logical_array = isnan(Theta);
[nan_indices] = find(nan_logical_array);
Theta(nan_indices) = 0;



u_vor = uazi.*sin(Theta) + U_field;
w_vor = -uazi.*cos(Theta);


figure
quiver((X(1:4:end,1:4:end))./r_omega,...
    (Z(1:4:end,1:4:end))./r_omega,...
    u_vor(1:4:end,1:4:end),w_vor(1:4:end,1:4:end)...
    ,2,'color',[1.00,0.00,0.00],...
    'MaxHeadSize',0.2,'LineWidth',1.5);
set(gca,'TickLabelInterpreter','latex','FontSize',13,...
    'XGrid','on','YGrid','on')
xlabel('x/r$_{\omega}$','Interpreter','Latex','FontSize',14);
ylabel('z/r$_{\omega}$','Interpreter','Latex','FontSize',14);
axis equal
xlim([-N-0.1 N+0.1])
ylim([-N-0.1 N+0.1])

Delx = X(1,2) - X(1,1);
Delz = Z(2,1) - Z(1,1);

[dudx, dudz] = gradient(u_vor,Delx,Delz);
[dwdx, dwdz] = gradient(w_vor,Delx,Delz);
omega = dwdx - dudz;
Mat = zeros(size(dudx,1),size(dudx,2),2,2);
Mat(:,:,1,1) = dudx;
Mat(:,:,1,2) = dudz;
Mat(:,:,2,1) = dwdx;
Mat(:,:,2,2) = dwdz;
Lambda_ci_Im_1 = zeros(size(u_vor));
Lambda_ci_Im_2 = zeros(size(u_vor));
Lambda_ci_Re_1 = zeros(size(u_vor));
Lambda_ci_Re_2 = zeros(size(u_vor));
for r = 1:size(Z,1)
    
    for c = 1:size(X,2)
        temp = reshape(Mat(r,c,:,:),[2,2]);
        eigval = eig(temp);
        Lambda_ci_Im_1(r,c) = imag(eigval(1,1));
        Lambda_ci_Im_2(r,c) = imag(eigval(2,1));
        Lambda_ci_Re_1(r,c) = real(eigval(1,1));
        Lambda_ci_Re_2(r,c) = real(eigval(2,1));
    end
end

figure
set(gcf,'Position',[764,373,906,686])
contourf(X/r_omega, Z/r_omega, Lambda_ci_Im_1, 501,'LineStyle','none')
cmap2 = getPyPlot_cMap('seismic', 51, false,...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
colormap(gca,cmap2)
caxis([-10 10])
hcb2=colorbar;
hcb2.Position=[0.889099628280219,0.109915542063659,0.022015823873409,0.801546391752577];
title(hcb2,'$Im(eig val(1))$','Interpreter','Latex','FontSize',15)
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation='out';
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('x/$r_{\omega}$','Interpreter','Latex','FontSize',14);
ylabel('z/$r_{\omega}$','Interpreter','Latex','FontSize',14);
axis equal
xlim([-N-0.1 N+0.1])
ylim([-N-0.1 N+0.1])

figure
set(gcf,'Position',[764,373,906,686])
contourf(X/r_omega, Z/r_omega, Lambda_ci_Im_2, 501,'LineStyle','none')
cmap2 = getPyPlot_cMap('seismic', 51, false,...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
colormap(gca,cmap2)
caxis([-10 10])
hcb2=colorbar;
hcb2.Position=[0.889099628280219,0.109915542063659,0.022015823873409,0.801546391752577];
title(hcb2,'$Im(eig val(2))$','Interpreter','Latex','FontSize',15)
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation='out';
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('x/$r_{\omega}$','Interpreter','Latex','FontSize',14);
ylabel('z/$r_{\omega}$','Interpreter','Latex','FontSize',14);
axis equal
xlim([-N-0.1 N+0.1])
ylim([-N-0.1 N+0.1])

figure
set(gcf,'Position',[764,373,906,686])
contourf(X/r_omega, Z/r_omega, Lambda_ci_Re_1, 501,'LineStyle','none')
cmap2 = getPyPlot_cMap('seismic', 51, false,...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
colormap(gca,cmap2)
caxis([-10 10])
hcb2=colorbar;
hcb2.Position=[0.889099628280219,0.109915542063659,0.022015823873409,0.801546391752577];
title(hcb2,'$Re(eig val(1))$','Interpreter','Latex','FontSize',15)
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation='out';
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('x/$r_{\omega}$','Interpreter','Latex','FontSize',14);
ylabel('z/$r_{\omega}$','Interpreter','Latex','FontSize',14);
axis equal
xlim([-N-0.1 N+0.1])
ylim([-N-0.1 N+0.1])

figure
set(gcf,'Position',[764,373,906,686])
contourf(X/r_omega, Z/r_omega, Lambda_ci_Re_2, 501,'LineStyle','none')
cmap2 = getPyPlot_cMap('seismic', 51, false,...
    '"C:\Users\ehsan010\Anaconda3\envs\Matlabpy\python.exe"');
colormap(gca,cmap2)
caxis([-10 10])
hcb2=colorbar;
hcb2.Position=[0.889099628280219,0.109915542063659,0.022015823873409,0.801546391752577];
title(hcb2,'$Re(eig val(2))$','Interpreter','Latex','FontSize',15)
hcb2.TickLabelInterpreter = 'latex';
hcb2.TickLength = 0;
hcb2.Label.Interpreter = 'latex';
hcb2.AxisLocation='out';
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('x/$r_{\omega}$','Interpreter','Latex','FontSize',14);
ylabel('z/$r_{\omega}$','Interpreter','Latex','FontSize',14);
axis equal
xlim([-N-0.1 N+0.1])
ylim([-N-0.1 N+0.1])