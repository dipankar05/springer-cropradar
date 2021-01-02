%% Modified Water Cloud Model: Semi-empirical model forms
% Prevot L, Champion I, Guyot G (1993) Estimating surface soil moisture and leaf area
%index of a wheat canopy using a dual-frequency (C and X bands) scatterometer.
%Remote Sens Environ 46(3):331–339
%%
clear;
%% Wheat: Freq: 9 GHz 
A = [0.021 0.417 0.99 19.2 -14.0]; 
thr = 40*pi/180;
%%-----------------------------------------

%%
%% Wheat: Freq: 5.3 GHz HH
% A = [0.0 0.089 0.0 19.2 -6.5]; % 
% thr = 20*pi/180;
%%
%% Alfalfa: Freq: 17 GHz 




%% Soil moisture (gravimetric) 
smg = 0.1; %0.08;

%% vegetation moisture fraction
% h = 2.5;

%%
%%
a = A(1); b = A(2); e = A(3); d = A(4); c = A(5); 
% % mp = 0.7;
lai = linspace(0.0,6.0,200);
% mv = 0.85.*lai; %%kg/m2
%[hg, smg] = meshgrid(h, sm);


%% Calibrated WCM
sigmasdb = c + d.*smg;
tau = exp((-2).*b.*lai.*sec(thr));
sigmac = a.*(lai.^e).*cos(thr).*(1-tau);
sigmas = (10.^((sigmasdb./10))).*tau;

sigma0g = sigmac + sigmas;

figure
x0=10;
y0=10;
width=700;
height=650;
set(gcf,'position',[x0,y0,width,height])
plot(lai,10.*log10(sigmac),'Color',[12 166 24]./255,'LineWidth',1.6)
hold on
% plot(lai,10.*log10(sigmast),'Color',[144, 252, 3]./255,'LineWidth',1.6)
plot(lai,10.*log10(sigmas),'Color',[148, 79, 4]./255,'LineWidth',1.6)
plot(lai,10.*log10(sigma0g),'black','LineWidth',2.0)
ylim([-35 0])
set(gca,'fontsize',24); 
legend('\sigma^0_{veg}','\sigma^0_{soil}','\sigma^0_{total}',...
    'fontsize',18,'interpreter','latex',...
    'Orientation','horizontal')
xlabel('LAI ($m^2~m^{-2}$)','fontsize',24,'interpreter','latex'); 
ylabel('$\sigma^0 (dB)$','fontsize',24,'interpreter','latex');

