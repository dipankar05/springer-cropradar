%%
clear;

%f = 1.8; % frequency in GHz
thetaid = 0: 0.5: 90; 
nf = length(thetaid);

%% Circular disk--Particle property
% % a = 2.25/100; % 2.25 cm radius
% % h = 0.5/1000; % 0.5 mm thickness
% % rho = 3000; % number density
% % d = 1.0; % canopy layer thickness 1.0m
% % er is calculated for mg = 0.5 g/g at frequency f GHz
%% Change the particle properties in the main function if required
%%

for ii = 1: nf
[sigma0vvdB(ii), sigma0hhdB(ii), sigma0vvdBr(ii), sigma0hhdBr(ii)] = Karam_Fung_GRG1988(thetaid(ii));
end

figure('DefaultAxesFontSize',14)
plot(thetaid', sigma0vvdB','r-','LineWidth',1.5)
grid on
hold on
plot(thetaid', sigma0hhdB','b-','LineWidth',1.5)
plot(thetaid', sigma0vvdBr','r--','LineWidth',1.5)
plot(thetaid', sigma0hhdBr','b--','LineWidth',1.5)
%plot(thetai', sigma0hvdB')
xlabel('\theta_i (deg.)')
ylabel('\sigma^\circ (dB)')
legend('\sigma^\circ_{VV}GRG','\sigma^\circ_{HH}GRG','\sigma^\circ_{VV}Rayleigh','\sigma^\circ_{HH}Rayleigh','Location','southwest','Orientation','vertical');
legend boxoff;
xlim([0 90]);
%ylim([-40 -15]);
% % Set graph size
x0=0;
y0=0;
width=350;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])

