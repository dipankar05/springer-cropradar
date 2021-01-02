%% S2RT model for layer of vegetation over a rough surface soil
% % Estimating backscatter cross sections
clc;
clear;

%f = 1.8; % frequency in GHz
% Radar incidence angle
thetai = 10: 1: 50; 
nf = length(thetai);

tic
%% Circular disk--Particle property
% % a = 3/100; % 3 cm radius
% % c = 0.1/1000; % 0.1 mm thickness
% % rho = 3000; % number density
% % h = 1.0; % canopy layer thickness 1.0m
% % er is calculated for mg = 0.5 g/g at frequency of 1.5 GHz
%% Change the particle properties in the main function if required
%%

for ii = 1: nf
[sigma0tvvdB(ii),sigma0gvvdB(ii),sigma0cvvdB(ii),sigma0cgvvdB(ii),sigma0gcgvvdB(ii),sigma0thhdB(ii),sigma0ghhdB(ii),sigma0chhdB(ii),sigma0cghhdB(ii),sigma0gcghhdB(ii)] = S2RT_Ulaby(thetai(ii));
fprintf('Incidence angle processing step: %d \n',ii);
end
t1 = cputime
%% Plotting co-polarized backscatter
%
figure('DefaultAxesFontSize',14)
plot(thetai', sigma0tvvdB','k-','LineWidth',1.8)
grid on
hold on
plot(thetai', sigma0gvvdB','Color' , [19 19 212]./255,'LineWidth',1.5)
plot(thetai', sigma0cvvdB','Color' , [21 181 10]./255,'LineWidth',1.5)
plot(thetai', sigma0cgvvdB','Color' , [240 10 10]./255,'LineWidth',1.5)
plot(thetai', sigma0gcgvvdB','Color' , [185 138 8]./255,'LineWidth',1.5)
%plot(thetai', sigma0hvdB')
xlabel('\theta_i (deg.)')
ylabel('\sigma^\circ_{VV} (dB)')
legend('\sigma^\circ_{t}','\sigma^\circ_{g}','\sigma^\circ_{c}','\sigma^\circ_{cg}','\sigma^\circ_{gcg}','Location','southwest','Orientation','vertical');
legend boxoff;
xlim([10 50]);
ylim([-45 -10]);
% % Set graph size
x0=0;
y0=0;
width=300;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])


figure('DefaultAxesFontSize',14)
plot(thetai', sigma0thhdB','k-','LineWidth',1.8)
grid on
hold on
plot(thetai', sigma0ghhdB','Color' , [19 19 212]./255,'LineWidth',1.5)
plot(thetai', sigma0chhdB','Color' , [21 181 10]./255,'LineWidth',1.5)
plot(thetai', sigma0cghhdB','Color' , [240 10 10]./255,'LineWidth',1.5)
plot(thetai', sigma0gcghhdB','Color' , [185 138 8]./255,'LineWidth',1.5)
%plot(thetai', sigma0hvdB')
xlabel('\theta_i (deg.)')
ylabel('\sigma^\circ_{HH} (dB)')
legend('\sigma^\circ_{t}','\sigma^\circ_{g}','\sigma^\circ_{c}','\sigma^\circ_{cg}','\sigma^\circ_{gcg}','Location','southwest','Orientation','vertical');
legend boxoff;
xlim([10 50]);
ylim([-45 -10]);
% % Set graph size
x0=0;
y0=0;
width=300;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])