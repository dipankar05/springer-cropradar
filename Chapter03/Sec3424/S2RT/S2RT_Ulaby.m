function [sigma0tvvdB,sigma0gvvdB,sigma0cvvdB,sigma0cgvvdB,sigma0gcgvvdB,sigma0thhdB,sigma0ghhdB,sigma0chhdB,sigma0cghhdB,sigma0gcghhdB] = S2RT_Ulaby(thetai)
%% 
f=1.5; % Frequency(GHz);
fg=f*1e9;
eo=single(8.854e-12); % permitivity of free space 
muo=single(pi*4e-7);  % permeability of free space
w=single(2*pi*fg);     % radial frequency
c=300e6;              % speed of light in free space
lambda0=c/fg;          % wave length in free space
ko=w*sqrt(muo*eo);    % wave number in free space
% er = 28.52 +2.13j;
%% Vegetation dielectric component ev (complex number)
% Particle moisture content mg (g/g)
mg = 0.7;
% Ulaby and El-Rayes (1987) Vegetation dielectric model based on Deybe-Cole dual-dispersion model
% [Ref] Ulaby, F.T. and El-Rayes, M.A., 1987. Microwave dielectric spectrum of vegetation-Part II: Dual-dispersion model. IEEE Transactions on Geoscience and Remote Sensing, (5), pp.550-557.
er = 1.7 - (0.74*mg) + (6.16*mg*mg);
vfw = mg*((0.55*mg)-0.076);
vb = (4.64*mg*mg)/(1+(7.36*mg*mg));
sigma = 1.27;
ev = er + (vfw*(4.9 + (75.0/(1+((1j*f)/18))) - ((1j*18*sigma)/f))) + (vb*(2.9 + (55.0/(1 + sqrt((1j*f)/0.18)))));
er = ev';
%%
%
% radar incidence angle
%thetai = 30;
theta0 = deg2rad(thetai);

%% Circular disk--Particle property
% % a = 10/1000; % 10 mm radius
% % c = 0.087/1000; % 0.087 mm thickness
% % V0 = (4/3)*pi*a*a*c;
% % n0 = 500; % number of particle
% % h = 1.0; % Canopy height in meter

a = 3.0/100; % 3.0 cm radius
c = 0.1/1000; % 0.1 mm thickness
V0 = (4/3)*pi*a*a*c;
n0 = 3000; % number of particle
h = 1.0; % Canopy height in meter

%% Calculation of a
Aa = pi/(2*a*a*a);
Ab = Aa;
Ac = 2/(a*a*c);
v = 0.5*a*a*c*(er-1);
at = (er-1)/(v*Aa +1);
an = (er-1)/(v*Ac +1);
C0 = ((w*w*muo*eo)/(4*pi))*V0;

%% Oriented volume
% % Distribution of gamma
lr = deg2rad(0);
ur = deg2rad(90);
dr = abs(ur-lr);
%% Scattering amplitude calculation
f2 = @(beta,gamma) (1/(pi/2)).*(1/dr).*(((0.5.*(sin(beta).^2).*(cos(theta0).^2) + (cos(beta).^2).*(sin(theta0).^2)).*(cos(gamma).^2)) + (0.5.*(cos(theta0).^2).*(sin(gamma).^2)));
f22 = integral2(f2,0,pi/2,lr,ur);
fvvii = C0*n0*(at + f22.*(an-at)); % <fvv(+i,+i)>
%--------------------------------------------
f3 = @(beta,gamma) (1/(pi/2)).*(1/dr).*((2.*(((0.5.*(sin(beta).^2).*(cos(theta0).^2) + (cos(beta).^2).*(sin(theta0).^2)).*(cos(gamma).^2)) + (0.5.*(cos(theta0).^2).*(sin(gamma).^2))).*real(at'.*(an-at))) + ((((3/8).*(cos(theta0).^4).*(sin(gamma).^4)) + (((sin(theta0).^4).*(cos(beta).^4) + (3/8).*(cos(theta0).^4).*(sin(beta).^4) + (3/16).*(sin(2.*theta0).^2).*(sin(2.*beta).^2)).*(cos(gamma).^4)) + ((3/16).*((cos(theta0).^4).*(sin(beta).^2) + (sin(2.*theta0).^2).*(cos(beta).^2)).*(sin(2.*gamma).^2))).* ((abs(an-at)).^2)));
f33 = integral2(f3,0,pi/2,lr,ur);
fvvii2 = C0*C0*n0*(((abs(at)).^2) + f33); % <|fvv(-i,i)|2>
%---------------------------------------------
f4 = @(beta,gamma) (1/(pi/2)).*(1/dr).*(0.5.*((sin(beta).^2).*(cos(gamma).^2) + (sin(gamma).^2)));
f44 = integral2(f4,0,pi/2,lr,ur);
fhhii = C0*n0*(at + f44.*(an-at)); % <fhh(i,i)>
%---------------------------------------------
f5 = @(beta,gamma) (1/(pi/2)).*(1/dr).*((((sin(beta).^2).*(cos(gamma).^2) + (sin(gamma).^2)).*real(at'.*(an-at))) + ((3/8).*((sin(beta).^4).*(cos(gamma).^4) + 0.5.*(sin(beta).^2).*(sin(2*gamma).^2) + (sin(gamma).^4)).*((abs(an-at)).^2)));
f55 = integral2(f5,0,pi/2,lr,ur);
fhhii2 = C0*C0*n0*(((abs(at)).^2) + f55); % <|fhh(-i,i)|2>
%---------------------------------------------
f6 = @(beta,gamma) (1/(pi/2)).*(1/dr).*((((cos(theta0).^2).*(sin(beta).^4) + (sin(2.*beta).^2).*(sin(theta0).^2)).*(cos(gamma).^4)) + ((cos(theta0).^2).*(sin(gamma).^4)) + ((0.5.*(sin(beta).^2).*(cos(theta0).^2) + (cos(beta).^2).*(sin(theta0).^2)).*(sin(2.*gamma).^2)));
f66 = integral2(f6,0,pi/2,lr,ur);
fvhii2 = C0*C0*(n0/8).*(f66.*((abs(an-at)).^2)); % <|fhh(-i,i)|2>

%---------------------------------------------
% sigma0vv = ((2.*pi)./lambda0).*fvvii2.*cos(theta0)./(2.*imag(fvvii));
% sigma0hh = ((2.*pi)./lambda0).*fhhii2.*cos(theta0)./(2.*imag(fhhii));
% sigma0vh = ((2.*pi)./lambda0).*fvhii2.*cos(theta0)./(imag(fhhii) + imag(fvvii));
% sigma0vvdB = 10.*log10(sigma0vv)
% sigma0hhdB = 10.*log10(sigma0hh)
% sigma0vhdB = 10.*log10(sigma0vh)
keh = imag(fvvii);
kev = imag(fhhii);
tauh = keh.*h./cos(theta0);
tauv = kev.*h./cos(theta0);
Upsilonh = exp((-1).*tauh);
Upsilonv = exp((-1).*tauv);

sigmavbackvv = ((2.*pi)./lambda0).*fvvii2;
sigmavbackhh = ((2.*pi)./lambda0).*fhhii2;
%----------------------------------------------
%
%
%
%%
%%
%% Bistatic scattering
thetaid = (pi./2)-theta0;
%% Scattering amplitude calculation
%--------------------------------------------
f7 = @(beta,gamma) (1/(pi/2)).*(1/dr).*((2.*(((0.5.*(sin(beta).^2).*(cos(thetaid).^2) + (cos(beta).^2).*(sin(thetaid).^2)).*(cos(gamma).^2)) + (0.5.*(cos(thetaid).^2).*(sin(gamma).^2))).*real(at'.*(an-at))) + ((((3/8).*(cos(thetaid).^4).*(sin(gamma).^4)) + (((sin(thetaid).^4).*(cos(beta).^4) + (3/8).*(cos(thetaid).^4).*(sin(beta).^4) + (3/16).*(sin(2.*thetaid).^2).*(sin(2.*beta).^2)).*(cos(gamma).^4)) + ((3/16).*((cos(thetaid).^4).*(sin(beta).^2) + (sin(2.*thetaid).^2).*(cos(beta).^2)).*(sin(2.*gamma).^2))).* ((abs(an-at)).^2)));
f77 = integral2(f7,0,pi/2,lr,ur);
fvvii3 = C0*C0*n0*(((abs(at)).^2) + f77); % <|fvv(-i,i)|2>

f8 = @(beta,gamma) (1/(pi/2)).*(1/dr).*((((sin(beta).^2).*(cos(gamma).^2) + (sin(gamma).^2)).*real(at'.*(an-at))) + ((3/8).*((sin(beta).^4).*(cos(gamma).^4) + 0.5.*(sin(beta).^2).*(sin(2*gamma).^2) + (sin(gamma).^4)).*((abs(an-at)).^2)));
f88 = integral2(f8,0,pi/2,lr,ur);
fhhii3 = C0*C0*n0*(((abs(at)).^2) + f88); % <|fhh(-i,i)|2>
%---------------------------------------------
sigmavbistvv = ((2.*pi)./lambda0).*fvvii3;
sigmavbisthh = ((2.*pi)./lambda0).*fhhii3;
%%
%
%
%% Backscattering Sigma calculation for direct canopy
sigma0cvv = (sigmavbackvv.*cos(theta0)./(2.*kev)).*(1 - Upsilonv.*Upsilonv);
sigma0chh = (sigmavbackhh.*cos(theta0)./(2.*keh)).*(1 - Upsilonh.*Upsilonh);
% sigma0vh = ((2.*pi)./lambda0).*fvhii2.*cos(theta0)./(imag(fhhii) + imag(fvvii));
sigma0cvvdB = 10.*log10(sigma0cvv);
sigma0chhdB = 10.*log10(sigma0chh);
% sigma0vhdB = 10.*log10(sigma0vh)

%% Backscattering Sigma calculation for ground-canopy and canopy-ground incoherent sum
%
%% Fresnel reflectivity calculation
% soil relative dielctric
mv0 = 0.16; %volumetric moisture
S = 0.306; C = 0.135; rho_b = 1.4; %loam soil parameters
[err, ei] = RelDielConst_Soil(f,27, rho_b, mv0, S, C); %get the dielectric
epsa = err -1i*ei;
%epsa = 6.28 - 1i * 1.53; % corresponding to mv = 0.16
% fresnel reflectivity
[gammav, gammah] = FresnelReflectivity(epsa,theta0);
%
%
sigma0cgvv = 2.*sigmavbistvv.*gammav.*Upsilonv.*Upsilonv.*h;
sigma0cghh = 2.*sigmavbisthh.*gammah.*Upsilonh.*Upsilonh.*h;
sigma0cgvvdB = 10.*log10(sigma0cgvv);
sigma0cghhdB = 10.*log10(sigma0cghh);

%% Backscattering Sigma calculation for ground-canopy-ground
sigma0gcgvv = gammav.*gammav.*Upsilonv.*Upsilonv.*sigma0cvv;
sigma0gcghh = gammah.*gammah.*Upsilonh.*Upsilonh.*sigma0chh;
sigma0gcgvvdB = 10.*log10(sigma0gcgvv);
sigma0gcghhdB = 10.*log10(sigma0gcghh);

%
%
%% Soil surface backscattering estimation with I2EM  model
sig1 = 0.5 / 100 ; % rms height in m
L1 = 5.0 / 100; % correl length in m
sp = 2; xx= 1.5; % selection of correl func
[sigma0gvvdB1,sigma0ghhdB1,sigma0ghvdB1] = I2EM_Backscatter_model(f, sig1, L1, theta0, epsa, sp, xx);
% linear scale conversion and multiplication with transmittivity
sigma0gvv = (10.^(sigma0gvvdB1./10)).*Upsilonv.*Upsilonv;
sigma0ghh = (10.^(sigma0ghhdB1./10)).*Upsilonh.*Upsilonh;
sigma0ghv = (10.^(sigma0ghvdB1./10)).*Upsilonh.*Upsilonv;
sigma0gvvdB = 10.*log10(sigma0gvv);
sigma0ghhdB = 10.*log10(sigma0ghh);
sigma0ghvdB = 10.*log10(sigma0ghv);

%% S2RT Total backscattering
sigma0tvv = sigma0gvv + sigma0cvv + sigma0cgvv + sigma0gcgvv;
sigma0thh = sigma0ghh + sigma0chh + sigma0cghh + sigma0gcghh;
sigma0tvvdB = 10.*log10(sigma0tvv);
sigma0thhdB = 10.*log10(sigma0thh);
end