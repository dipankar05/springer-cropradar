%%
clear;

thetai = 0 : 0.5 : 85; 
nf = length(thetai);

%% Circular disk--Particle property

%%
for ii = 1: nf
[sigma0vvdB(ii), sigma0hhdB(ii)] = LeVine1984_PO_circulardisc(thetai(ii));
end

figure
plot(thetai', sigma0vvdB','r-','LineWidth',1.5)
grid on
hold on
plot(thetai', sigma0hhdB','b-','LineWidth',1.5)
%plot(thetai', sigma0hvdB')
xlabel('\theta_i (deg.)')
ylabel('\sigma^\circ (dB)')
legend('\sigma^\circ_{VV}','\sigma^\circ_{HH}','Location','southwest','Orientation','vertical');
legend boxoff;
% xlim([0 90]);
% ylim([-30 -10]);
% % Set graph size
% x0=0;
% y0=0;
% width=200;
% height=500;
% set(gcf,'units','points','position',[x0,y0,width,height])