%%
clear;

f = 1.8; % frequency in GHz
thetai = 0: 5: 85; 
nf = length(thetai);

%% Circular disk--Particle property
% % a = 2.25/100; % 2.25 cm radius
% % h = 0.5/1000; % 0.5 mm thickness
% % rho = 3000; % number density
% % d = 1.0; % canopy layer thickness 1.0m
% % er is calculated for mg = 0.5 g/g at frequency f GHz
%% Change the particle properties in the main function if required
%%
%% Orientation angle lower and upper range in degree 0<theta<180
lr = 60;
ur = 90;

%%
for ii = 1: nf
[sigma0vvdB(ii), sigma0hhdB(ii), sigma0hvdB(ii)] = Foldy_distortedBornApproximation_orienteddisk_LangSidhu1983(f,thetai(ii),lr,ur);
end

figure('DefaultAxesFontSize',14)
plot(thetai', sigma0vvdB','r-','LineWidth',1.5)
grid on
hold on
plot(thetai', sigma0hhdB','b-','LineWidth',1.5)
%plot(thetai', sigma0hvdB')
xlabel('\theta_i (deg.)')
ylabel('\sigma^\circ (dB)')
legend('\sigma^\circ_{VV}','\sigma^\circ_{HH}','Location','southwest','Orientation','vertical');
legend boxoff;
xlim([0 90]);
ylim([-30 -10]);
% % Set graph size
x0=0;
y0=0;
width=250;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])
