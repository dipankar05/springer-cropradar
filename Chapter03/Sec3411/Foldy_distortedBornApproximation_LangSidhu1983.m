%% Backscattering from a layer of vegetation with a discrete approach
% % Wave theory (Foldy-distorted Born Approximation)
% % Vegetation layer--> circular discs with uniform orientation

%% Ref: Land and Sidhu, "Electromagnetic Backscattering From a Layer of Vegetation: A Discrete Approach", 
% %...IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. GE-21, NO. 1, JANUARY 1983.

%%
function [sigma0vvdB sigma0hhdB sigma0hvdB] = FdBLangSidhu(f,thetai)
% Frequency (GHz)
%f = 1.8;
% radar incidence angle
%thetai = 25; %degree
theta0 = deg2rad(thetai);

%% Vegetation dielectric component ev (complex number)
% Particle moisture content mg (g/g)
mg = 0.5;
% Ulaby and El-Rayes (1987) Vegetation dielectric model based on Deybe-Cole dual-dispersion model
% [Ref] Ulaby, F.T. and El-Rayes, M.A., 1987. Microwave dielectric spectrum of vegetation-Part II: Dual-dispersion model. IEEE Transactions on Geoscience and Remote Sensing, (5), pp.550-557.
er = 1.7 - (0.74*mg) + (6.16*mg*mg);
vfw = mg*((0.55*mg)-0.076);
vb = (4.64*mg*mg)/(1+(7.36*mg*mg));
sigma = 1.27;
ev = er + (vfw*(4.9 + (75.0/(1+((1j*f)/18))) - ((1j*18*sigma)/f))) + (vb*(2.9 + (55.0/(1 + sqrt((1j*f)/0.18)))));
er = ev';


%------------
C = 3*(10e+8); % velocity of light m/sec
lambda = C/(f*(10e+9));
k0 = (2*pi)/lambda;  % wavenumber
%---------------------
D = er -1;
ar = D/(1+D);
at = D;
ap = D;
%% Circular disk--Particle property
a = 2.25/100; % 2.25 cm radius
h = 0.5/1000; % 0.5 mm thickness
Vp = pi*a*a*h;
rho = 3000; % number density
delta = rho * Vp;
d = 1.0; % canopy layer thickness 1.0m

%% scattering amplitude calculation
% Particle orientation distribution theta -->0 to pi and phi -->0 to 2*pi

%% angular average of a parameters
% calculation axy
f1 = @(th,ph) (1/pi).*(1/(2*pi)).*(((ar.*(sin(th).^2) + at.*(cos(th).^2)) - ap).*cos(ph).*sin(ph));
axy = integral2(f1,0,pi,0,2*pi);
% calculation axz
f2 = @(th,ph) (1/pi).*(1/(2*pi)).*((ar-at).*sin(th).*cos(th).*cos(ph));
axz = integral2(f2,0,pi,0,2*pi);
% calculation ayz
f3 = @(th,ph) (1/pi).*(1/(2*pi)).*((ar-at).*sin(th).*cos(th).*sin(ph));
ayz = integral2(f3,0,pi,0,2*pi);
% calculation axx and ayy, azz
f4 = @(th,ph) (1/pi).*(1/(2*pi)).*(((ar.*(sin(th).^2) + at.*(cos(th).^2)).*(cos(ph).^2)) + (ap.*(sin(ph).^2)));
axx = integral2(f4,0,pi,0,2*pi);
f5 = @(th,ph) (1/pi).*(1/(2*pi)).*(((ar.*(sin(th).^2) + at.*(cos(th).^2)).*(sin(ph).^2)) + (ap.*(cos(ph).^2)));
ayy = integral2(f5,0,pi,0,2*pi);
f6 = @(th,ph) (1/pi).*(1/(2*pi)).*(ar.*(cos(th).^2) + at.*(sin(th).^2));
azz = integral2(f6,0,pi,0,2*pi);

%% Calculation of polarization terms with respect to a
ahh2 = (abs(ayy)).^2;
ahv2 = (cos(theta0).^2 .* (abs(axy)).^2) + (sin(theta0).^2 .* (abs(ayz)).^2);
avv2 = (cos(theta0).^4 .* (abs(axx)).^2) + (sin(theta0).^4 .* (abs(azz)).^2) + (cos(theta0).^2 .* sin(theta0).^2 .* ((4.*(abs(axz)).^2)+(2.* real(axx.*azz'))));  

%% calculation of extinction coefficients
kz0 = k0 .* cos(theta0);
kzh = kz0 + ((delta.*k0.*k0.*ayy)./(2.*kz0));
kzv = kz0 + (((delta.*k0.*k0)./(2.*kz0)).* (axx.*(cos(theta0).^2) + azz.*(sin(theta0).^2)));

%% Backscatter coefficient calculation
sigma0vv = (((k0.^4).*delta.*Vp)./(4*pi)).*avv2.*((1 - exp((-4).*imag(kzv).*d))./(4.*imag(kzv)));
sigma0hh = (((k0.^4).*delta.*Vp)./(4*pi)).*ahh2.*((1 - exp((-4).*imag(kzh).*d))./(4.*imag(kzh)));
sigma0hv = (((k0.^4).*delta.*Vp)./(4*pi)).*ahv2.*((1 - exp((-2).*imag(kzh+kzv).*d))./(2.*imag(kzh+kzv)));

sigma0vvdB = 10.*log10(sigma0vv);
sigma0hhdB = 10.*log10(sigma0hh);
sigma0hvdB = 10.*log10(sigma0hv);

end