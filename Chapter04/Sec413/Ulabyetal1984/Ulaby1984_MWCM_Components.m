%% Modified Water Cloud Model: Semi-empirical model forms
% Ulaby F, Allen C, Eger Iii G, Kanemasu E (1984) Relating the microwave backscattering
% coefficient to leaf area index. Remote Sens Environ 14(1-3):113–133
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
A = [0.218 2.56 0.411 0.025 0.0 0.197]; % @30 degree inc
thr = 50*pi/180;
%%
%% Alfalfa: Freq: 17 GHz 




%% Soil moisture (gravimetric) 
smg = 0.1; %0.08;

%% vegetation moisture fraction
h = 2.5;

%%
%%
al = A(1); bl = A(2); alpl = A(3); ast = A(4); alpst = A(5); cs = A(6);
% % mp = 0.7;
lai = linspace(0.0,2,200);
mv = 0.85.*lai; %%kg/m2
%[hg, smg] = meshgrid(h, sm);


%% Calibrated WCM

sigmal = al.*(1-exp((-bl).*lai./h)).*(1-exp((-2).*alpl.*sec(thr).*lai)).*cos(thr);
sigmast = ast.*mv.*h.*exp((-2).*alpl.*sec(thr).*lai);
sigmas = cs.*smg.*exp((-2).*alpl.*sec(thr).*lai).*exp((-2).*alpst.*mv.*h);
sigma0g = sigmal + sigmast + sigmas; 

figure
x0=10;
y0=10;
width=700;
height=650;
set(gcf,'position',[x0,y0,width,height])
plot(lai,10.*log10(sigmal),'Color',[12 166 24]./255,'LineWidth',1.6)
hold on
plot(lai,10.*log10(sigmast),'Color',[144, 252, 3]./255,'LineWidth',1.6)
plot(lai,10.*log10(sigmas),'Color',[148, 79, 4]./255,'LineWidth',1.6)
plot(lai,10.*log10(sigma0g),'black','LineWidth',2.0)
% ylim([-35 5])
set(gca,'fontsize',24); 
legend('\sigma^0_{l}','\sigma^0_{st}','\sigma^0_{soil}','\sigma^0_{total}',...
    'fontsize',18,'interpreter','latex',...
    'Orientation','horizontal')
xlabel('LAI ($m^2~m^{-2}$)','fontsize',24,'interpreter','latex'); 
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

