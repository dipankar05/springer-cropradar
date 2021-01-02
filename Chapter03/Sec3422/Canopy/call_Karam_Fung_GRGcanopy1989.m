%%
clear;

%f = 1.8; % frequency in GHz
f = 1.0: 0.1: 3; 
%theta0 = 25: 5: 55; 
nf = length(f);

%% Circular disk--Particle property
% % a = 2.25/100; % 2.25 cm radius
% % h = 0.5/1000; % 0.5 mm thickness
% % rho = 3000; % number density
% % d = 1.0; % canopy layer thickness 1.0m
% % er is calculated for mg = 0.5 g/g at frequency f GHz
%% Change the particle properties in the main function if required
%%

for ii = 1: nf
[sigma0vvdB(ii), sigma0hhdB(ii), sigma0vhdB(ii), sigma0hvdB(ii)] = Karam_Fung_GRGcanopy1989(f(ii));
end

figure('DefaultAxesFontSize',14)
% figure
plot(f', sigma0vvdB','r-','LineWidth',1.5)
grid on
hold on
plot(f', sigma0hhdB','b-','LineWidth',1.5)
plot(f', sigma0vhdB','g-','LineWidth',1.5)
%plot(thetai', sigma0hvdB')
xlabel('Frequency (GHz)')
ylabel('\sigma^\circ (dB)')
legend('\sigma^\circ_{VV}','\sigma^\circ_{HH}','\sigma^\circ_{VH}','Location','southwest','Orientation','vertical');
legend boxoff;
xlim([1 3]);
ylim([-50 -10]);
% % Set graph size
x0=0;
y0=0;
width=250;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])

