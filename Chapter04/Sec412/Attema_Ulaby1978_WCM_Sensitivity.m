%% Water Cloud Model: Semi-empirical model forms
% Attema E, Ulaby FT (1978) Vegetation modeled as a water cloud. Radio Sci 13(2):357–364
%%
clear;
%% Alfalfa: Freq: 8.6 GHz 
% A = [0.023 0.228 1.08]; % @0 degree inc
% thr = 0*pi/180;
%%-----------------------------------------
% A = [0.00603 0.228 1.08]; % @30 degree inc
% thr = 30*pi/180;
%%
%% Alfalfa: Freq: 13.0 GHz 
% A = [0.045 0.311 2.54]; % @0 degree inc
% thr = 0*pi/180;
%%-----------------------------------------
A = [0.009 0.311 1.54]; % @30 degree inc
thr = 30*pi/180;
%%
%%------------------------------------------------------
%%
a = A(1); c = A(2); d = A(3);  
% % mp = 0.7;
w = linspace(0,4.0,200);
h = 0.6;
sm = linspace(20,400,300);
[wg, smg] = meshgrid(w, sm);

%% Calibrated WCM

sigmac = c.*(1-(1.*exp((-d).*wg.*h.*sec(thr)))).*cos(thr);
sigmas = a.*smg.*cos(thr).*exp((-d).*wg.*h.*sec(thr));
sigma0g = sigmas + sigmac; 
% contour(hg, smg, sigma0g,1.2:-0.0001:0);
figure
x0=10;
y0=10;
width=700;
height=650;
set(gcf,'position',[x0,y0,width,height])
s3=surf(wg, smg, 10.*log10(sigma0g),'FaceAlpha',1);
colormap(jet)
colorbar
caxis([-10 10])
s3.EdgeColor = 'none';
set(gca,'fontsize',24); 
colorbar('FontSize',24,'FontName','Arial');
xlabel('Vegetation water content ($kg~m^{-3}$)','fontsize',24,'interpreter','latex'); 
ylabel('Volumetric soil moisture ($kg~m^{-3}$)','fontsize',24,'interpreter','latex');
axis square; box on;
view(2)

