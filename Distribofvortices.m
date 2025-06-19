Data2 = load('VF_PIVWT7.mat');
VF_PIVWT7 = Data2.VF_PIVWT7;

%% Profile of Lambda_{T}
WT7VorX_scales = load('G:\My Drive\Research\DATABASE\Vortex statistics\data\mesh_07ms_scales.mat');


figure
plot(WT7VorX_scales.lambda(1:22)*100,WT7VorX_scales.z(1:22)/WT7VorX_scales.delta)
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('$\lambda_{T}$[cm]','Interpreter','Latex','FontSize',14);
ylabel('$z/\delta$','Interpreter','Latex','FontSize',14);

Lambda_T_bar = mean(WT7VorX_scales.lambda(1:22),1);
std(WT7VorX_scales.lambda(1:22),0,1)

% %% Distribution of vortices in Lambda_{T} square
% 
% ix = 1;
% DistribofVorX = zeros(1,2);
% for S = 1:size(VF_PIVWT7.Lambda_ci,3)
%     
%     xint = floor(VF_PIVWT7.x(end)/Lambda_T_bar);
%     zint = floor(VF_PIVWT7.z(end)/Lambda_T_bar);
%     for r = 1 : zint
%         for c = 1 : xint
%             Matrix = VF_PIVWT7.Lambda_ci((r-1)*zint+1:(r)*zint+1,...
%                 (c-1)*xint+1:(c)*xint+1,S);            
%             binaryMatrix = Matrix ~= 0;
%             [labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);
%             
%             Numb_pro = 
%             Numb_ret = w
%             DistribofVorX (ix, :) = [Numb_pro, Numb_ret];
%             ix = ix + 1;
%             end
%         end
%     end
% end
