%%
clear;

%f = 1.8; % frequency in GHz
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

for ii = 1: nf
[sigma0vvdB(ii), sigma0hhdB(ii), sigma0vhdB(ii)] = Karam_Fung_RTM_CircularDiskOR_1983(thetai(ii));
end

figure('DefaultAxesFontSize',14)
plot(thetai', sigma0vvdB','r-','LineWidth',1.5)
grid on
hold on
plot(thetai', sigma0hhdB','b-','LineWidth',1.5)
plot(thetai', sigma0vhdB','g-','LineWidth',1.5)
%plot(thetai', sigma0hvdB')
xlabel('\theta_i (deg.)')
ylabel('\sigma^\circ (dB)')
legend('\sigma^\circ_{VV}','\sigma^\circ_{HH}','\sigma^\circ_{VH}','Location','southwest','Orientation','vertical');
legend boxoff;
xlim([0 90]);
ylim([-40 -15]);
% % Set graph size
x0=0;
y0=0;
width=250;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])

