function [sigma0vvdB sigma0hhdB] = POLeVine1984(thetai)
%% Physical optics model
f = 9.0; %GHz
%thetai = 20;
theta0 = deg2rad(thetai);
a = 20/100;
T = 5/100;

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
%er = 11+5j;

%------------
C = 3*(10e+8); % velocity of light m/sec
lambda = C/(f*(10e+9));
k0 = (2*pi)/lambda;  % wavenumber
%------------
I = 0 +1j;
G = sqrt(er - ((sin(theta0)).^2));
zitah = 1;
zitav = 1/er;
rh = (cos(theta0) - zitah.*G)./(cos(theta0) + zitah.*G);
rv = (cos(theta0) - zitav.*G)./(cos(theta0) + zitav.*G);
th = (2.*cos(theta0).*sqrt(zitah))./(cos(theta0) + zitah.*G);
tv = (2.*cos(theta0).*sqrt(zitav))./(cos(theta0) + zitav.*G);
alpha = 1./(2.*k0.*T);
S0 = 1./sqrt(4.*pi).*T.*k0.*k0.*(er-1);
Op = cos(theta0) + G; % Omega +ve
Om = cos(theta0) - G; % Omega -ve
gp = ((sin(theta0)).^2) - cos(theta0).*G; % gamma +ve
gm = ((sin(theta0)).^2) + cos(theta0).*G; % gamma -ve
ev = tv.*exp((-1).*I.*alpha.*Om)./(1 - rv.*rv.*exp(I.*4.*alpha.*G));
eh = th.*exp((-1).*I.*alpha.*Om)./(1 - rh.*rh.*exp(I.*4.*alpha.*G));
Q = 2.*k0.*a.*sin(theta0);
Svt = pi.*a.*a.*besselj(1,Q)./(0.5.*Q);

sigma0hh = ((abs((sinc(alpha.*Op) - (rh.*exp(I.*2.*alpha.*G).*sinc(alpha.*Om))).*eh.*S0)).^2).*Svt.*Svt;
sigma0vv = ((abs((gp.*sinc(alpha.*Op) - (rv.*exp(I.*2.*alpha.*G).*gm.*sinc(alpha.*Om))).*(ev./sqrt(er)).*S0)).^2).*Svt.*Svt;

sigma0hhdB = 10.*log10(sigma0hh);
sigma0vvdB = 10.*log10(sigma0vv);
end

