%% Dielectric-Slab model: Semi-empirical model forms
% Bush T, Ulaby F (1976) Radar return from a continuous vegetation canopy. IEEE Trans Antennas Propag 24(3):269–276
%%
clear;
%% Alfalfa: Freq: 8.6 GHz 
% A = [0.469 9.941 26.73 0.35];
%% Alfalfa: Freq: 13.0 GHz 
% A = [0.269 12.827 46.01 0.46];
%% Alfalfa: Freq: 17.0 GHz 
A = [0.556 6.585 35.01 0.47];
%%------------------------------------------------------

%% Soil moisture (gravimetric) 
smg = 0.08;

%% vegetation moisture fraction
mp = 0.7;
%%
a = A(1); b = A(2); c = A(3); d = A(4); 
% plant height
hg = linspace(0.0,0.8,200);

%[hg, smg] = meshgrid(h, sm);

%% Calibrated WCM

sigmac = (d.*sqrt(mp).*hg);
sigmas = a.*exp((b.*smg)- (c.*sqrt(mp).*(hg.^2.6)));
sigma0g = sigmas + sigmac; 

figure
x0=10;
y0=10;
width=700;
height=650;
set(gcf,'position',[x0,y0,width,height])
plot(hg,10.*log10(sigmac),'Color',[12 166 24]./255,'LineWidth',1.6)
hold on
plot(hg,10.*log10(sigmas),'Color',[148, 79, 4]./255,'LineWidth',1.6)
plot(hg,10.*log10(sigma0g),'black','LineWidth',2.0)
ylim([-35 5])
set(gca,'fontsize',24); 
legend('\sigma^0_{veg}','\sigma^0_{soil}','\sigma^0_{total}',...
    'fontsize',18,'interpreter','latex',...
    'Orientation','horizontal')
xlabel('$h~(m)$','fontsize',24,'interpreter','latex'); 
ylabel('$\sigma^0 (dB)$','fontsize',24,'interpreter','latex');



