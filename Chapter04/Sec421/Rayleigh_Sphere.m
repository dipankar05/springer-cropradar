%% Rayleigh scattering --Spherical particle
clear;
clc;
f=5.0; % Frequency(GHz);
fq=f*1e9;
eo=single(8.854e-12); % permitivity of free space 
muo=single(pi*4e-7);  % permeability of free space
w=single(2*pi*fq);     % radial frequency
c0=300e6;              % speed of light in free space
lambda0=c0/fq;          % wave length in free space
k0=w*sqrt(muo*eo);    % wave number in free space
%er = 35.8 + 3.54i;
%% Vegetation dielectric component ev (complex number)
% Particle moisture content mg (g/g)
mg = 0.2:0.1:0.9;
n0 = 200:400:4200;
[x,y] = size(mg);
[xx,yy] = size(n0);
for ik = 1:y
% Ulaby and El-Rayes (1987) Vegetation dielectric model based on Deybe-Cole dual-dispersion model
% [Ref] Ulaby, F.T. and El-Rayes, M.A., 1987. Microwave dielectric spectrum of vegetation-Part II: Dual-dispersion model. IEEE Transactions on Geoscience and Remote Sensing, (5), pp.550-557.
err = 1.7 - (0.74.*mg(ik)) + (6.16.*mg(ik).*mg(ik));
vfw = mg(ik).*((0.55.*mg(ik))-0.076);
vb = (4.64.*mg(ik).*mg(ik))/(1+(7.36.*mg(ik).*mg(ik)));
sigma = 1.27;
ev(ik) = err + (vfw*(4.9 + (75.0/(1+((1j*f)/18))) - ((1j*18*sigma)/f))) + (vb*(2.9 + (55.0/(1 + sqrt((1j*f)/0.18)))));
er(ik) = ev(ik);
% Vegetation water content VWC (kg/m3)
rhow = 500; % particle wet density, kg/m3
vwc1(ik) = mg(ik)*rhow;
%%
%% Circular disk--Particle property
r0 = 1.0/1000; % radius of sphere r mm
v0 = (4/3).*pi.*r0.^3;
h = 2.0; % vegetation height
%% Scattering and extinction amplitudes
S = 128.*(pi.^5).*(r0.^6).*(abs((1-er(ik))./(2+er(ik))).^2)./(3.*(lambda0.^4));
k = 8.*(pi.^2).*(r0.^3).*imag((1-er(ik))./(2+er(ik)))./lambda0;

for ii=1:yy
%
%% Final backscatter crosssection
ke(ik,ii) = n0(ii).*k;
sigmavbackhh(ik,ii) = n0(ii).*S;
%------------------------------------------------------------------------------------
% VWC sphere:
vwc(ik,ii) = vwc1(ik)*v0*n0(ii)*h;
%------------------------------------------------------------------------------------
% WCM parameters
Bgrghh(ik,ii) = 2.*ke(ik,ii).*h;
Agrghh(ik,ii) = sigmavbackhh(ik,ii)./2.*ke(ik,ii);

end
end



%%
%%
%
%
for ik = 1:y
% Ulaby and El-Rayes (1987) Vegetation dielectric model based on Deybe-Cole dual-dispersion model
% [Ref] Ulaby, F.T. and El-Rayes, M.A., 1987. Microwave dielectric spectrum of vegetation-Part II: Dual-dispersion model. IEEE Transactions on Geoscience and Remote Sensing, (5), pp.550-557.
err = 1.7 - (0.74.*mg(ik)) + (6.16.*mg(ik).*mg(ik));
vfw = mg(ik).*((0.55.*mg(ik))-0.076);
vb = (4.64.*mg(ik).*mg(ik))/(1+(7.36.*mg(ik).*mg(ik)));
sigma = 1.27;
ev(ik) = err + (vfw*(4.9 + (75.0/(1+((1j*f)/18))) - ((1j*18*sigma)/f))) + (vb*(2.9 + (55.0/(1 + sqrt((1j*f)/0.18)))));
er(ik) = ev(ik);
% Vegetation water content VWC (kg/m3)
rhow = 500; % particle wet density, kg/m3
vwc1(ik) = mg(ik)*rhow;
%%
%% Circular disk--Particle property
r0 = 3.0/1000; % radius of sphere r mm
v0 = (4/3).*pi.*r0.^3;
h = 2.0; % vegetation height
%% Scattering and extinction amplitudes
S = 128.*(pi.^5).*(r0.^6).*(abs((1-er(ik))./(2+er(ik))).^2)./(3.*(lambda0.^4));
k = 8.*(pi.^2).*(r0.^3).*imag((1-er(ik))./(2+er(ik)))./lambda0;

for ii=1:yy
%
%% Final backscatter crosssection
ke(ik,ii) = n0(ii).*k;
sigmavbackhh(ik,ii) = n0(ii).*S;
%------------------------------------------------------------------------------------
% VWC sphere:
vwc3(ik,ii) = vwc1(ik)*v0*n0(ii)*h;
%------------------------------------------------------------------------------------
% WCM parameters
Bgrghh3(ik,ii) = 2.*ke(ik,ii).*h;
Agrghh3(ik,ii) = sigmavbackhh(ik,ii)./2.*ke(ik,ii);

end
end

%%
%
%
%
%%
%%
%
%
for ik = 1:y
% Ulaby and El-Rayes (1987) Vegetation dielectric model based on Deybe-Cole dual-dispersion model
% [Ref] Ulaby, F.T. and El-Rayes, M.A., 1987. Microwave dielectric spectrum of vegetation-Part II: Dual-dispersion model. IEEE Transactions on Geoscience and Remote Sensing, (5), pp.550-557.
err = 1.7 - (0.74.*mg(ik)) + (6.16.*mg(ik).*mg(ik));
vfw = mg(ik).*((0.55.*mg(ik))-0.076);
vb = (4.64.*mg(ik).*mg(ik))/(1+(7.36.*mg(ik).*mg(ik)));
sigma = 1.27;
ev(ik) = err + (vfw*(4.9 + (75.0/(1+((1j*f)/18))) - ((1j*18*sigma)/f))) + (vb*(2.9 + (55.0/(1 + sqrt((1j*f)/0.18)))));
er(ik) = ev(ik);
% Vegetation water content VWC (kg/m3)
rhow = 500; % particle wet density, kg/m3
vwc1(ik) = mg(ik)*rhow;
%%
%% Circular disk--Particle property
r0 = 5.0/1000; % radius of sphere r mm
v0 = (4/3).*pi.*(r0.^3);
h = 2.0; % vegetation height
%% Scattering and extinction amplitudes
S = 128.*(pi.^5).*(r0.^6).*((abs((1-er(ik))./(2+er(ik)))).^2)./(3.*(lambda0.^4));
k = 8.*(pi.^2).*(r0.^3).*imag((1-er(ik))./(2+er(ik)))./lambda0;

for ii=1:yy
%
%% Final backscatter crosssection
ke(ik,ii) = n0(ii).*k;
sigmavbackhh(ik,ii) = n0(ii).*S;
%------------------------------------------------------------------------------------
% VWC sphere:
vwc5(ik,ii) = vwc1(ik)*v0*n0(ii)*h;
%------------------------------------------------------------------------------------
% WCM parameters
Bgrghh5(ik,ii) = 2.*ke(ik,ii).*h;
Agrghh5(ik,ii) = sigmavbackhh(ik,ii)./2.*ke(ik,ii);

end
end



%% Plotting mg vs WCM parameter
% Regression fitting:
dd2 = mg';
gg2 = Agrghh(:,3);
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit2 = (abs(coeff(1)*dd2.^coeff(2)));
% %
gg2 = Agrghh3(:,3);
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit3 = (abs(coeff(1)*dd2.^coeff(2)));
% %
gg2 = Agrghh5(:,3);
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit5 = (abs(coeff(1)*dd2.^coeff(2)));

figure('DefaultAxesFontSize',14)
for kk=1:yy
% Circular disk
%---------------------------------------------------------------------------------
% create a default color map ranging from blue to light blue or cyan
length = yy;
blue = [0, 0, 1];
cyan = [149, 235, 255]./255;
% RGB Table help: https://www.rapidtables.com/web/color/RGB_Color.html
colors_p = [linspace(blue(1),cyan(1),length)', linspace(blue(2),cyan(2),length)', linspace(blue(3),cyan(3),length)'];
%----------------
green = [0,144,10]./255;
g1 = [0,1,0];
colors_g = [linspace(green(1),g1(1),length)', linspace(green(2),g1(2),length)', linspace(green(3),g1(3),length)'];
%----------------
red = [168,7,7]./255;
r1 = [255,94,0]./255;
colors_r = [linspace(red(1),r1(1),length)', linspace(red(2),r1(2),length)', linspace(red(3),r1(3),length)'];
%----------------
loglog(mg,abs(Agrghh(:,kk)'),'p','MarkerEdgeColor','None','MarkerFaceColor',colors_p(kk,:));
hold on
loglog(mg,abs(Agrghh3(:,kk)'),'o','MarkerEdgeColor',colors_g(kk,:),'MarkerFaceColor',colors_g(kk,:));
loglog(mg,abs(Agrghh5(:,kk)'),'d','MarkerEdgeColor','None','MarkerFaceColor',colors_r(kk,:));
grid on
%---------------------------------------------------------------------------------
loglog(mg',Bfit2,'k-','LineWidth',0.5)
loglog(mg',Bfit3,'k-','LineWidth',0.5)
loglog(mg',Bfit5,'k-','LineWidth',0.5)
%---------------------------------------------------------------------------------
xlabel('Particle moisture content $m_g (g~g^{-1})$','Interpreter','latex');
ylabel('$\sigma_{v,HH}/2k_{e,HH}$','Interpreter','latex');
xlim([0 1.0]);
end
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])

clear Bfit2; 
clear Bfit3; 
clear Bfit5; 
%% Ploting VWC vs WCM parameter
% Regression fitting:
dd = log10(vwc(8,:)');
gg = log10(Agrghh(8,:)');
P3 = polyfit(dd,gg,1);  % y = p1*x + p2
Bfit = 10.^(P3(1).*dd+P3(2));
%
dd = log10(vwc3(8,:)');
gg = log10(abs(Agrghh3(8,:))');
P3 = polyfit(dd,gg,1); % y = p1*x + p2
Bfit3 = 10.^(P3(1)*dd+P3(2));
%
dd = log10(vwc5(8,:)');
gg = log10(abs(Agrghh5(8,:))');
P3 = polyfit(dd,gg,1); % y = p1*x + p2
Bfit5 = 10.^(P3(1)*dd+P3(2));
%
%
figure('DefaultAxesFontSize',14)
loglog(vwc(:),abs(Agrghh(:)),'p','MarkerEdgeColor','None','MarkerFaceColor',[0 0 1])
hold on
loglog(vwc3(:),abs(Agrghh3(:)),'o','MarkerEdgeColor','None','MarkerFaceColor',[0 1 0])
loglog(vwc5(:),abs(Agrghh5(:)),'d','MarkerEdgeColor','None','MarkerFaceColor',[1 0 0])
xlabel('Vegetation water content VWC, $(kg~m^{-3})$','Interpreter','latex');
ylabel('$\sigma_{v,HH}/2k_{e,HH}$','Interpreter','latex');
loglog(vwc(8,:)',Bfit,'k-','LineWidth',0.5)
loglog(vwc3(8,:)',Bfit3,'k-','LineWidth',0.5)
loglog(vwc5(8,:)',Bfit5,'k-','LineWidth',0.5)
grid on
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])

clear Bfit; 
clear Bfit2; 
clear Bfit3; 
clear Bfit5; 
%% Plotting B parameter
%% Plotting mg vs WCM parameter
% Regression fitting:
dd2 = mg';
gg2 = Bgrghh(:,3);
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit2 = (abs(coeff(1)*dd2.^coeff(2)));
%
gg2 = Bgrghh3(:,3);
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit3 = (abs(coeff(1)*dd2.^coeff(2)));
%
gg2 = Bgrghh5(:,3);
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit5 = (abs(coeff(1)*dd2.^coeff(2)));
%
%
figure('DefaultAxesFontSize',14)
for kk=1:yy
% Circular disk
%---------------------------------------------------------------------------------
% create a default color map ranging from blue to light blue or cyan
length = yy;
blue = [0, 0, 1];
cyan = [149, 235, 255]/255;
% RGB Table help: https://www.rapidtables.com/web/color/RGB_Color.html
colors_p = [linspace(blue(1),cyan(1),length)', linspace(blue(2),cyan(2),length)', linspace(blue(3),cyan(3),length)'];
loglog(mg,abs(Bgrghh(:,kk)'),'p','MarkerEdgeColor','None','MarkerFaceColor',colors_p(kk,:))
hold on
loglog(mg,abs(Bgrghh3(:,kk)'),'o','MarkerEdgeColor','None','MarkerFaceColor',colors_g(kk,:))
loglog(mg,abs(Bgrghh5(:,kk)'),'d','MarkerEdgeColor','None','MarkerFaceColor',colors_r(kk,:))
grid on
xlabel('Particle moisture content $m_g (g~g^{-1})$','Interpreter','latex');
ylabel('$k_{e,HH}.h$','Interpreter','latex');
%xlim([0 1.0]);
%---------------------------------------------------------------------------------
loglog(mg',Bfit2,'k-','LineWidth',0.5)
loglog(mg',Bfit3,'k-','LineWidth',0.5)
loglog(mg',Bfit5,'k-','LineWidth',0.5)
end
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])

clear Bfit;
clear Bfit2; 
clear Bfit3; 
clear Bfit5; 
%% Ploting VWC vs WCM parameter
% Regression fitting:
dd = log10(vwc(8,:)');
gg = log10(abs(Bgrghh(8,:)'));
P3 = polyfit(dd,gg,1);  % y = p1*x + p2
Bfit = 10.^(P3(1)*dd+P3(2));
%
dd = log10(vwc3(8,:)');
gg = log10(abs(Bgrghh3(8,:)'));
P3 = polyfit(dd,gg,1);  % y = p1*x + p2
Bfit3 = 10.^(P3(1)*dd+P3(2));
%
dd = log10(vwc5(8,:)');
gg = log10(abs(Bgrghh5(8,:)'));
P3 = polyfit(dd,gg,1);  % y = p1*x + p2
Bfit5 = 10.^(P3(1)*dd+P3(2));
%
figure('DefaultAxesFontSize',14)
loglog(vwc(:),abs(Bgrghh(:)),'p','MarkerEdgeColor','None','MarkerFaceColor',[0 0 1])
hold on
loglog(vwc3(:),abs(Bgrghh3(:)),'o','MarkerEdgeColor','None','MarkerFaceColor',[0 1 0])
loglog(vwc5(:),abs(Bgrghh5(:)),'d','MarkerEdgeColor','None','MarkerFaceColor',[1 0 0])
grid on
% legend('Circular disk','Needle/Cylinder','Location','northwest','Orientation','vertical');
% legend boxoff;
xlabel('Vegetation water content VWC, $(kg~m^{-3})$','Interpreter','latex');
ylabel('$k_{e,HH}.h$','Interpreter','latex');
loglog(vwc(8,:)',Bfit,'k-','LineWidth',0.5)
loglog(vwc3(8,:)',Bfit3,'k-','LineWidth',0.5)
loglog(vwc5(8,:)',Bfit5,'k-','LineWidth',0.5)
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])
%