%% Dielectric-Slab model: Semi-empirical model forms
% Bush T, Ulaby F (1976) Radar return from a continuous vegetation canopy. IEEE Trans Antennas Propag 24(3):269–276
%%
clear;
%% Alfalfa: Freq: 8.6 GHz 
A = [0.469 9.941 26.73 0.35]; % Obtained by fitting data
%% Alfalfa: Freq: 13.0 GHz 
% A = [0.269 12.827 46.01 0.46];
%% Alfalfa: Freq: 17.0 GHz 
% A = [0.556 6.585 35.01 0.47];
%%------------------------------------------------------
%%
a = A(1); b = A(2); c = A(3); d = A(4); 
mp = 0.7;
h = linspace(0.0,0.8,200);
sm = linspace(0.05,0.3,200);
[hg, smg] = meshgrid(h, sm);

%% Calibrated WCM

sigmac = (d.*sqrt(mp).*hg);
sigmas = a.*exp((b.*smg)- (c.*sqrt(mp).*(hg.^2.6)));
sigma0g = sigmas + sigmac; 
% contour(hg, smg, sigma0g,1.2:-0.0001:0);
figure
x0=10;
y0=10;
width=700;
height=650;
set(gcf,'position',[x0,y0,width,height])
s3=surf(hg, smg, 10.*log10(sigma0g),'FaceAlpha',1);
colormap(jet)
colorbar
%caxis([-2 10])
s3.EdgeColor = 'none';
set(gca,'fontsize',24); 
colorbar('FontSize',24,'FontName','Arial');
xlabel('$h (m)$','fontsize',24,'interpreter','latex'); 
ylabel('$m_s (gm~cm^{-3})$','fontsize',24,'interpreter','latex');
axis square; box on;
view(2)

