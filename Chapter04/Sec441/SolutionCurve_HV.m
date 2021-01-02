% Visualization of climatological and radiometric bounds
clc; clear; close all;

%Radar incident angle=35
thr=(35*pi)/180;
%Inputs, tbl=[LAI,SM]
A = [0.03857004,  0.55096747,  0.02209724, -0.01739473,  0.33006099];
a = A(1); b = A(2); c = A(3); d = A(4); e = A(5);
% fHH = ((a*x1.^e).*((1-exp((-2).*b.*x1))*sec(thr))) + ...
%     ((c + d*x2).*exp((-2).*b.*x1*sec(thr)));

lai = linspace(0.01,0.65,200);
sm = linspace(0.05,0.50,200);
[laig, smg] = meshgrid(lai, sm);

%% Calibrated WCM

fHHg = ((a*laig.^e).*((1-exp((-2).*b.*laig))*sec(thr))) + (c +(d*smg).*exp((-2).*b.*laig*sec(thr)));

% cost function
func_H = (fHHg - 0.025).^2;

% func_H = (fHHg - 0.02).^2;


%% plot
hh=figure('Color',[1 1 1]);
% hh = figure
x0=10;
y0=10;
width=700;
height=650;
set(gcf,'position',[x0,y0,width,height]);
cmap = [204,236,230;153,216,201;102,194,164;65,174,118;35,139,69;0,109,44;0,68,27]./255;
colormap(cmap);

% 121
hold on;
% contour(laig, smg, func_H,0.0001:-0.00000001:0);
contour(laig, smg, func_H,0.005:-0.000001:0);
set(gca,'fontsize',16); 
idx = find(func_H<=0.000000005);
pt = plot(laig(idx), smg(idx),'b-','LineWidth',1.5);
colorbar('FontSize',24,'FontName','Arial');
%%
% plot([max(laig(idx)) max(laig(idx))],[0.05 0.5],'--','lineWidth',1.5,'Color',[67,67,67]./255);
% plot([0 0.65],[min(smg(idx)) min(smg(idx))],'--','lineWidth',1.5,'Color',[67,67,67]./255);
%%
xlabel('LAI $(m^2~ m^{-2})$','fontsize',18,'interpreter','latex'); 
ylabel('m$_v (m^3~ m^{-3})$','fontsize',18,'interpreter','latex');
axis square; box on; 
% xlim([0 6.5]); 
ylim([0.05 0.50]); h.YTick = [0.1 0.2 0.3 0.4 0.5];
legend([pt],{'Solution Curve'},...
     'Orientation','horizontal','Box','off','FontSize',16);
 set(gca,'fontsize',24);