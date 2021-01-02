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
%% Alfalfa: Freq: 17 GHz 




%% Soil moisture (gravimetric) 
smg = 80; %0.08;

%% vegetation moisture fraction
h = 0.6;
%%
%%
a = A(1); c = A(2); d = A(3);  
% % mp = 0.7;
wg = linspace(0.0,4.0,200);

%[hg, smg] = meshgrid(h, sm);


%% Calibrated WCM

sigmac = c.*(1-(exp((-d).*wg.*h/sec(thr))))*cos(thr);
sigmas = a.*smg*cos(thr).*exp((-d).*wg.*h.*sec(thr));
sigma0g = sigmas + sigmac; 

figure
x0=10;
y0=10;
width=700;
height=650;
set(gcf,'position',[x0,y0,width,height])
plot(wg,10.*log10(sigmac),'Color',[12 166 24]./255,'LineWidth',1.6)
hold on
plot(wg,10.*log10(sigmas),'Color',[148, 79, 4]./255,'LineWidth',1.6)
plot(wg,10.*log10(sigma0g),'black','LineWidth',2.0)
ylim([-35 5])
set(gca,'fontsize',24); 
legend('\sigma^0_{veg}','\sigma^0_{soil}','\sigma^0_{total}',...
    'fontsize',18,'interpreter','latex',...
    'Orientation','horizontal')
xlabel('Vegetation water content ($kg~m^{-3}$)','fontsize',24,'interpreter','latex'); 
ylabel('$\sigma^0 (dB)$','fontsize',24,'interpreter','latex');

% % contour(hg, smg, sigma0g,1.2:-0.0001:0);
% s3=surf(hg, smg, sigma0g,'FaceAlpha',1);
% colormap(jet)
% colorbar
% caxis([-2 10])
% s3.EdgeColor = 'none';

% colorbar('FontSize',10,'FontName','Arial');
% xlabel('$h (m)$','fontsize',18,'interpreter','latex'); 
% ylabel('$m_s (gm~cm^{-3})$','fontsize',18,'interpreter','latex');
% axis square; box on;
% view(2)

