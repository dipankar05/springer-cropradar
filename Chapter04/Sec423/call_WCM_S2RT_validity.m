%% S2RT model for layer of vegetation over a rough surface soil
% % Estimating backscatter cross sections
clc;
clear;

%f = 1.8; % frequency in GHz
% Radar incidence angle
thetai = 15: 5: 45; 
nf = length(thetai);
% Particle moisture content mg (g/g)
mg = 0.1:0.05:0.9;
yy = length(mg);
%tic
%% Circular disk--Particle property
% % a = 3/100; % 3 cm radius
% % c = 0.1/1000; % 0.1 mm thickness
% % rho = 3000; % number density
% % h = 1.0; % canopy layer thickness 1.0m
% % er is calculated for mg = 0.5 g/g at frequency of 1.5 GHz
%% Change the particle properties in the main function if required
%%
for jj = 1:yy
for ii = 1: nf
[sigma0tvvedB(jj,ii),sigma0thhedB(jj,ii)] = S2RT_Ulaby(mg(jj),thetai(ii));
fprintf('Incidence angle processing step: %d \n',ii);
end
fprintf('Mg processing step: %d \n',jj);
end
%t1 = cputime;


%% Plotting co-polarized backscatter
% VV plot
data = sigma0tvvedB;
a = mean(data,2);
mina = min(data,[],2);
maxa = max(data,[],2);
A = [mg',mina,a,maxa];

figure('DefaultAxesFontSize',14)
plot(A(:,1), A(:,3), 'sr', 'MarkerFaceColor','r')    
hold on
plot([A(:,1) A(:,1)]', [A(:,2) A(:,4)]','-r')   
hold on
plot(mg',a,'-r','LineWidth',0.5)
xlabel('Particle moisture content $m_g (g~g^{-1})$','Interpreter','latex');
ylabel('\sigma^\circ_{total}- \sigma^\circ_{WCM}(dB)')
% legend('\sigma^\circ_{t}','\sigma^\circ_{g}','\sigma^\circ_{c}','\sigma^\circ_{cg}','\sigma^\circ_{gcg}','Location','southwest','Orientation','vertical');
% legend boxoff;
% % Set graph size
x0=0;
y0=0;
width=500;
height=200;
set(gcf,'units','points','position',[x0,y0,width,height])

%% HH Plot
clear data;
clear a;
clear A;
data = sigma0thhedB;
a = mean(data,2);
mina = min(data,[],2);
maxa = max(data,[],2);
A = [mg',mina,a,maxa];

figure('DefaultAxesFontSize',14)
plot(A(:,1), A(:,3), 'sb', 'MarkerFaceColor','b')    
hold on
plot([A(:,1) A(:,1)]', [A(:,2) A(:,4)]','-b')   
hold on
plot(mg',a,'-b','LineWidth',0.5)
xlabel('Particle moisture content $m_g (g~g^{-1})$','Interpreter','latex');
ylabel('\sigma^\circ_{total}- \sigma^\circ_{WCM}(dB)')
% legend('\sigma^\circ_{t}','\sigma^\circ_{g}','\sigma^\circ_{c}','\sigma^\circ_{cg}','\sigma^\circ_{gcg}','Location','southwest','Orientation','vertical');
% legend boxoff;
ylim([0 0.3]);
% % Set graph size
x0=0;
y0=0;
width=500;
height=200;
set(gcf,'units','points','position',[x0,y0,width,height])

