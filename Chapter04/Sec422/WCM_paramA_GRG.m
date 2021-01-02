%%
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
mg = 0.1:0.1:0.9;
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
er(ik) = ev(ik)';
% Vegetation water content VWC (kg/m3)
rhow = 500; % particle wet density, kg/m3
vwc1(ik) = mg(ik)*rhow;
%%
%% Circular disk--Particle property
% Fig-8 Karam 1989
% a = 10/1000; % 10 mm radius
a = 3.0/100; % 2.5 cm radius
t = 0.1/1000; % 0.5 mm thickness
v0 = (4/3)*pi*a*a*t;
%n0 = 500; % number of particle density
h = 2.0; % Canopy height in meter

m = 2*a/t;
% Calculating demagnetizing factors and modifying functions
% [Ref] Karam, M.A., Fung, A.K. and Antar, Y.M., 1988. Electromagnetic wave scattering from some vegetation samples. IEEE Transactions on Geoscience and Remote Sensing, 26(6), pp.799-808.

gt = (0.5/((m*m)-1))*((((m*m)/sqrt((m*m)-1))*asin(sqrt((m*m)-1)/m))-1);
gn = ((m*m)/((m*m)-1))*(1-((1/sqrt((m*m)-1))*asin(sqrt((m*m)-1)/m)));
at = 1/(((er(ik)-1)*gt)+1);
an = 1/(((er(ik)-1)*gn)+1);
%%
% costhetasd = ((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis));
% 
% costhetaid = (cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii));
% 
% cosphisd = (((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (costhetasd.^2)));
%
% cosphiid = (((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - (costhetaid.^2)));
%
% sinphisd = (sqrt(1 - (cosphisd.^2)));
% sinphiid = (sqrt(1 - (cosphiid.^2)));
% sinthetaid = (sqrt(1 - (costhetaid.^2)));
% sinthetasd = (sqrt(1 - (costhetasd.^2)));
%% Debye modifying function
% % qxa = (k0.*(sin(thetai).*cos(alpha-phii) - sin(thetas).*cos(alpha-phis)));
% % qya = (k0.*(sin(thetai).*sin(alpha-phii) - sin(thetas).*sin(alpha-phis)));
% % qza = (k0.*(cos(thetai)+cos(thetas)));
% % Qe = (a.*sqrt(((cos(beta).*qxa + sin(beta).*qza).^2) + (qya.^2)));
% % Qe = (a.*sqrt(((cos(beta).*(k0.*(sin(thetai).*cos(alpha-phii) - sin(thetas).*cos(alpha-phis))) + sin(beta).*(k0.*(cos(thetai)+cos(thetas)))).^2) + ((k0.*(sin(thetai).*sin(alpha-phii) - sin(thetas).*sin(alpha-phis))).^2)));
% % % mus = (2 * besselj(1,Qe))./Qe;
% % mus = ((2 * besselj(1,(a.*sqrt(((cos(beta).*(k0.*(sin(thetai).*cos(alpha-phii) - sin(thetas).*cos(alpha-phis))) + sin(beta).*(k0.*(cos(thetai)+cos(thetas)))).^2) + ((k0.*(sin(thetai).*sin(alpha-phii) - sin(thetas).*sin(alpha-phis))).^2)))))./(a.*sqrt(((cos(beta).*(k0.*(sin(thetai).*cos(alpha-phii) - sin(thetas).*cos(alpha-phis))) + sin(beta).*(k0.*(cos(thetai)+cos(thetas)))).^2) + ((k0.*(sin(thetai).*sin(alpha-phii) - sin(thetas).*sin(alpha-phis))).^2))));
%%
% %% Backscattering geometry
gamma = 0;
theta0 = 35;
thetai = deg2rad(theta0);
thetas = thetai;
phii = 0;
phis = phii+pi;
%% Range of beta and probability
lr = deg2rad(0);
ur = deg2rad(30);
Pb = ur-lr;


%% Scattering amplitude calculation
% % mus-----------integral
f7 = @(alpha,beta) (1).*(((a.*sqrt(((cos(beta).*(k0.*(sin(thetai).*cos(alpha-phii) - sin(thetas).*cos(alpha-phis))) + sin(beta).*(k0.*(cos(thetai)+cos(thetas)))).^2) + ((k0.*(sin(thetai).*sin(alpha-phii) - sin(thetas).*sin(alpha-phis))).^2)))));
%f77 = integral2(f7,0,2*pi,0,deg2rad(30));
mus = @(alpha,beta) 2.*besselj(1,f7(alpha,beta))./f7(alpha,beta); % GRG Approximation
% mus = @(alpha,beta) 1; % In case of Rayleigh scattering
%%
% % HH-------------------------------
%f2 = @(alpha,beta) (1/(2*pi)).*(1/(pi)).*(1).*(cosphisd.*cosphiid - sinphisd.*sinphiid).*at.*mus;
f2 = @(alpha,beta) (1).*(((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))) .*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2)))) - (sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2)) .*sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2)))).*at.*mus(alpha,beta);
% f22 = integral2(f2,0,2*pi,0,pi);
fhhsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f2(alpha,beta); % <fhh(-i,i)>


% % VV-------------------------
%f3 = @(alpha,beta) (1/(2*pi)).*(1/(pi)).*(1).*((an.*sinthetaid.*sinthetasd) - (at.*costhetaid.*costhetasd.*(cosphisd.*cosphiid + sinphisd.*sinphiid))).*mus;
f3 = @(alpha,beta) (1).*((an.*(sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).*(sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2)))) - (at.*(cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).*((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).*((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))) + (sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2))).*(sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2)))))).*(mus(alpha,beta));
%f33 = integral2(f3,0,2*pi,0,deg2rad(30));
fvvsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f3(alpha,beta); % <fhh(-i,i)>

% % VH---------------------------
%f6 = @(alpha,beta) (1/(2*pi)).*(1/(pi/2)).*(1).*(costhetasd.*(sinphisd.*cosphiid - cosphisd.*sinphiid)).*at.*mus;
f6 = @(alpha,beta) (1).*(((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).*((sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2))).*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))) - (((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).*(sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2))))).*at.*mus(alpha,beta);
%f66 = integral2(f6,0,2*pi,0,deg2rad(30));
fvhsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f6(alpha,beta); % <fvh(-i,i)>

% % HV--------------------------
%f7 = @(alpha,beta) (1/(2*pi)).*(1/(pi/2)).*(1).*(costhetaid.*(sinphisd.*cosphiid - cosphisd.*sinphiid)).*at.*mus;
f7 = @(alpha,beta) (1).*((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).*((sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2))).*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))) - (((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).*(sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2))))).*at.*mus(alpha,beta);
fhvsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f7(alpha,beta); % <fhv(-i,i)>

%%
%%f10 = @(alpha,beta) f8(alpha,beta).*f9(alpha,beta);
tvi = @(alpha,beta) (-1).*sin(beta).*cos(thetai).*cos(alpha-phii) - cos(beta).*sin(thetai);
thi = @(alpha,beta) sin(beta).*sin(alpha-phii);
tvs = @(alpha,beta) sin(beta).*cos(thetas).*cos(alpha-phis) + cos(beta).*sin(thetas);
ths = @(alpha,beta) sin(beta).*sin(alpha-phis);
Dsi = @(alpha,beta) sqrt((((tvs(alpha,beta)).^2) + ((ths(alpha,beta)).^2)).*(((tvi(alpha,beta)).^2) + ((thi(alpha,beta)).^2)));


%% Backscatter calculation in reference frame
% % VV-----------------------
F10 = @(alpha,beta) (1./Dsi(alpha,beta)).*((tvs(alpha,beta).* (fvvsi(alpha,beta).*tvi(alpha,beta) - fvhsi(alpha,beta).*thi(alpha,beta))) - (ths(alpha,beta).* (fhvsi(alpha,beta).*tvi(alpha,beta) - fhhsi(alpha,beta).*thi(alpha,beta))));
F101 = (1/(2*pi)).*(1./Pb).*(1).*(integral2(F10,0,2*pi,lr,ur));

% % HV-----------------------
F11 = @(alpha,beta) (1./Dsi(alpha,beta)).*((ths(alpha,beta).* (fvvsi(alpha,beta).*tvi(alpha,beta) - fvhsi(alpha,beta).*thi(alpha,beta))) + (tvs(alpha,beta).* (fhvsi(alpha,beta).*tvi(alpha,beta) - fhhsi(alpha,beta).*thi(alpha,beta))));
F111 = (1/(2*pi)).*(1./Pb).*(1).*(integral2(F11,0,2*pi,lr,ur));

% % VH-----------------------
F12 = @(alpha,beta) (1./Dsi(alpha,beta)).*((tvs(alpha,beta).* (fvvsi(alpha,beta).*thi(alpha,beta) + fvhsi(alpha,beta).*tvi(alpha,beta))) - (ths(alpha,beta).* (fhvsi(alpha,beta).*thi(alpha,beta) + fhhsi(alpha,beta).*tvi(alpha,beta))));
F121 = (1/(2*pi)).*(1./Pb).*(1).*(integral2(F12,0,2*pi,lr,ur));

% % HH-----------------------
F13 = @(alpha,beta) (1./Dsi(alpha,beta)).*((ths(alpha,beta).* (fvvsi(alpha,beta).*thi(alpha,beta) + fvhsi(alpha,beta).*tvi(alpha,beta))) + (tvs(alpha,beta).* (fhvsi(alpha,beta).*thi(alpha,beta) + fhhsi(alpha,beta).*tvi(alpha,beta))));
F131 = (1/(2*pi)).*(1./Pb).*(1).*(integral2(F13,0,2*pi,lr,ur));

%%

%%
% 
%% Forward scattering geometry
%gamma = 0;
%thetai = deg2rad(30);
thetas = pi - thetai;
%phii = 0;
phis = phii;
% 
% In case of forward scattering Debye function mus =1.0;
% 
%% Scattering amplitude calculation

% % mus-----------integral

mus = @(alpha,beta) 1;
%%
% % HH-------------------------------
%f2 = @(alpha,beta) (1/(2*pi)).*(1/(pi)).*(1).*(cosphisd.*cosphiid - sinphisd.*sinphiid).*at.*mus;
f2 = @(alpha,beta) (1).*(((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))) .*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2)))) - (sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2)) .*sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2)))).*at.*mus(alpha,beta);
% f22 = integral2(f2,0,2*pi,0,pi);
fhhsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f2(alpha,beta); % <fhh(-i,i)>


% % VV-------------------------
%f3 = @(alpha,beta) (1/(2*pi)).*(1/(pi)).*(1).*((an.*sinthetaid.*sinthetasd) - (at.*costhetaid.*costhetasd.*(cosphisd.*cosphiid + sinphisd.*sinphiid))).*mus;
f3 = @(alpha,beta) (1).*((an.*(sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).*(sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2)))) - (at.*(cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).*((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).*((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))) + (sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2))).*(sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2)))))).*(mus(alpha,beta));
%f33 = integral2(f3,0,2*pi,0,deg2rad(30));
fvvsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f3(alpha,beta); % <fhh(-i,i)>

% % VH---------------------------
%f6 = @(alpha,beta) (1/(2*pi)).*(1/(pi/2)).*(1).*(costhetasd.*(sinphisd.*cosphiid - cosphisd.*sinphiid)).*at.*mus;
f6 = @(alpha,beta) (1).*(((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).*((sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2))).*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))) - (((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).*(sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2))))).*at.*mus(alpha,beta);
%f66 = integral2(f6,0,2*pi,0,deg2rad(30));
fvhsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f6(alpha,beta); % <fvh(-i,i)>

% % HV--------------------------
%f7 = @(alpha,beta) (1/(2*pi)).*(1/(pi/2)).*(1).*(costhetaid.*(sinphisd.*cosphiid - cosphisd.*sinphiid)).*at.*mus;
f7 = @(alpha,beta) (1).*((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).*((sqrt(1 - ((((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).^2))).*(((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))) - (((sin(thetas).*(cos(gamma).*cos(beta).*cos(alpha-phis) - sin(gamma).*sin(alpha-phis))) - (cos(thetas).*cos(gamma).*sin(beta)))./sqrt(1 - (((-1).*cos(beta).*cos(thetas) - sin(beta).*sin(thetas).*cos(alpha-phis)).^2))).*(sqrt(1 - ((((sin(thetai).*(cos(gamma).*cos(beta).*cos(alpha-phii) - sin(gamma).*sin(alpha-phii))) + (cos(thetai).*cos(gamma).*sin(beta)))./sqrt(1 - ((cos(beta).*cos(thetai) - sin(beta).*sin(thetai).*cos(alpha-phii)).^2))).^2))))).*at.*mus(alpha,beta);
fhvsi = @(alpha,beta) ((k0.*k0.*(er(ik)-1).*v0)./(4*pi)).*f7(alpha,beta); % <fhv(-i,i)>

%%
%%f10 = @(alpha,beta) f8(alpha,beta).*f9(alpha,beta);
tvi = @(alpha,beta) (-1).*sin(beta).*cos(thetai).*cos(alpha-phii) - cos(beta).*sin(thetai);
thi = @(alpha,beta) sin(beta).*sin(alpha-phii);
tvs = @(alpha,beta) sin(beta).*cos(thetas).*cos(alpha-phis) + cos(beta).*sin(thetas);
ths = @(alpha,beta) sin(beta).*sin(alpha-phis);
Dsi = @(alpha,beta) sqrt((((tvs(alpha,beta)).^2) + ((ths(alpha,beta)).^2)).*(((tvi(alpha,beta)).^2) + ((thi(alpha,beta)).^2)));


%% Backscatter calculation in reference frame
% % VV-----------------------
F14 = @(alpha,beta) (1./Dsi(alpha,beta)).*((tvs(alpha,beta).* (fvvsi(alpha,beta).*tvi(alpha,beta) - fvhsi(alpha,beta).*thi(alpha,beta))) - (ths(alpha,beta).* (fhvsi(alpha,beta).*tvi(alpha,beta) - fhhsi(alpha,beta).*thi(alpha,beta))));
F141 = (1/(2*pi)).*(1./Pb).*(1).*(integral2(F14,0,2*pi,lr,ur));

% % HH-----------------------
F15 = @(alpha,beta) (1./Dsi(alpha,beta)).*((ths(alpha,beta).* (fvvsi(alpha,beta).*thi(alpha,beta) + fvhsi(alpha,beta).*tvi(alpha,beta))) + (tvs(alpha,beta).* (fhvsi(alpha,beta).*thi(alpha,beta) + fhhsi(alpha,beta).*tvi(alpha,beta))));
F151 = (1/(2*pi)).*(1./Pb).*(1).*(integral2(F15,0,2*pi,lr,ur));

%%
% 
for ii=1:yy
%
%% Final backscatter crosssection
ksv = n0(ii).*n0(ii).*2.*2.*pi.*((abs(F101).^2) + (abs(F111).^2));
ksh = n0(ii).*n0(ii).*2.*2.*pi.*((abs(F121).^2) + (abs(F131).^2));
% HH--------------------------------------
% FHHF for forward scattering case
% FHH for backward scattering case
FHHF = n0(ii).*F151;
FHH = n0(ii).*F131;
kah = n0(ii).*abs((4*pi/k0)*imag(FHHF));
khh(ik,ii) = ksh + kah;
sigmavbackhh(ik,ii) = 4.*pi.*(abs(FHH).^2);

% VV--------------------------------------
FVVF = n0(ii).*F141;
FVV = n0(ii).*F101;
kav = n0(ii).*abs((4*pi/k0)*imag(FVVF));
kvv(ik,ii) = ksv + kav;
sigmavbackvv(ik,ii) = 4.*pi.*(abs(FVV).^2);

%------------------------------------------------------------------------------------
% VWC circular disk:
vwc(ik,ii) = vwc1(ik)*v0*n0(ii)*h;
%------------------------------------------------------------------------------------
% WCM parameters
Bgrgvv(ik,ii) = 2.*kvv(ik,ii).*h;
Agrgvv(ik,ii) = sigmavbackvv(ik,ii)./2.*kvv(ik,ii);
Bgrghh(ik,ii) = 2.*khh(ik,ii).*h;
Agrghh(ik,ii) = sigmavbackhh(ik,ii)./2.*khh(ik,ii); 

end
end

%% Plotting
% 
% Regression fitting:
dd2 = mg';
gg2 = abs(Agrgvv(:,1));
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit2 = log10(abs(coeff(1)*dd2.^coeff(2)));
%Bfit2 = log10(abs(2e-7*dd2.^1.2869));

%% Plotting mg vs WCM parameter
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
loglog(mg,abs(Agrgvv(:,kk)'),'o','MarkerEdgeColor','None','MarkerFaceColor',colors_p(kk,:))
hold on
grid on
xlabel('Particle moisture content $m_g (g~g^{-1})$','Interpreter','latex');
ylabel('$\sigma_{v,VV}/2k_{e,VV}$','Interpreter','latex');
%xlim([0 1.0]);
loglog(mg',10.^Bfit2,'k-','LineWidth',0.5)
end
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])


%% Ploting VWC vs WCM parameter
% Regression fitting:
dd = (vwc(:));
gg = (abs(Agrgvv(:)));
% P = polyfit(dd,gg,1)  % y = p1*x + p2
% Bfit = (P(1).*dd+P(2));
f1 = fit(dd,gg,'power1');
coeff=coeffvalues(f1);
Bfit = (abs(coeff(1)*dd.^coeff(2)));

figure('DefaultAxesFontSize',14)
loglog(vwc(:),abs(Agrgvv(:)),'o','MarkerEdgeColor','None','MarkerFaceColor',[0 0 1])
hold on
grid on
% legend('Circular disk','Needle/Cylinder','Location','northwest','Orientation','vertical');
% legend boxoff;
xlabel('Vegetation water content VWC, $(kg~m^{-3})$','Interpreter','latex');
ylabel('$\sigma_{v,VV}/2k_{e,VV}$','Interpreter','latex');
loglog(vwc(:),Bfit,'k-','LineWidth',0.5)
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])

%% HH polarization
%
%
%% Plotting mg vs WCM parameter
% Regression fitting:
dd2 = mg';
gg2 = Agrghh(:,1);
f1 = fit(dd2,gg2,'power1');
coeff=coeffvalues(f1);
Bfit2 = log10(abs(coeff(1)*dd2.^coeff(2)));
%Bfit2 = log10(abs(2e-7*dd2.^1.2869));

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
loglog(mg,abs(Agrghh(:,kk)'),'o','MarkerEdgeColor','None','MarkerFaceColor',colors_p(kk,:))
hold on
grid on
xlabel('Particle moisture content $m_g (g~g^{-1})$','Interpreter','latex');
ylabel('$\sigma_{v,HH}/2k_{e,HH}$','Interpreter','latex');
%xlim([0 1.0]);
%---------------------------------------------------------------------------------
loglog(mg',10.^Bfit2,'k-','LineWidth',0.5)
end
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])


%% Ploting VWC vs WCM parameter
% Regression fitting:
dd = (vwc(:));
gg = (abs(Agrghh(:)));
% P = polyfit(dd,gg,1)  % y = p1*x + p2
% Bfit = (P(1).*dd+P(2));
f1 = fit(dd,gg,'power1');
coeff=coeffvalues(f1);
Bfit = (abs(coeff(1)*dd.^coeff(2)));

figure('DefaultAxesFontSize',14)
loglog(vwc(:),abs(Agrghh(:)),'o','MarkerEdgeColor','None','MarkerFaceColor',[0 0 1])
hold on
grid on
xlabel('Vegetation water content VWC, $(kg~m^{-3})$','Interpreter','latex');
ylabel('$\sigma_{v,HH}/2k_{e,HH}$','Interpreter','latex');
loglog(vwc(:),Bfit,'k-','LineWidth',0.5)
% % Set graph size
x0=0;
y0=0;
width=400;
height=350;
set(gcf,'units','points','position',[x0,y0,width,height])