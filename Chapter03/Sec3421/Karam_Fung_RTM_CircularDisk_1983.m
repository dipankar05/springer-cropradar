function [sigma0vvdB, sigma0hhdB, sigma0vhdB] = KaramFungRTMCD(thetai)
f=1.5; % Frequency(GHz);
f=f*1e9;
eo=single(8.854e-12); % permitivity of free space 
muo=single(pi*4e-7);  % permeability of free space
w=single(2*pi*f);     % radial frequency
c=300e6;              % speed of light in free space
lambda0=c/f;          % wave length in free space
ko=w*sqrt(muo*eo);    % wave number in free space
er = 28.52 +2.13j;
% radar incidence angle
theta0 = deg2rad(thetai);

%% Circular disk--Particle property
a = 10/1000; % 10 mm radius
c = 0.087/1000; % 0.087 mm thickness
V0 = (4/3)*pi*a*a*c;
n0 = 500; % number of particle
%% Calculation of a
Aa = pi/(2*a*a*a);
Ab = Aa;
Ac = 2/(a*a*c);
v = 0.5*a*a*c*(er-1);
at = (er-1)/(v*Aa +1);
an = (er-1)/(v*Ac +1);
C0 = ((w*w*muo*eo)/(4*pi))*V0;

%% Scattering amplitude calculation
f2 = @(beta,gamma) (1/(pi/2)).*(1/(pi/2)).*(((0.5.*(sin(beta).^2).*(cos(theta0).^2) + (cos(beta).^2).*(sin(theta0).^2)).*(cos(gamma).^2)) + (0.5.*(cos(theta0).^2).*(sin(gamma).^2)));
f22 = integral2(f2,0,pi/2,0,pi/2);
fvvii = C0*n0*(at + f22.*(an-at)); % <fvv(+-i,+-i)>
%--------------------------------------------
f3 = @(beta,gamma) (1/(pi/2)).*(1/(pi/2)).*((2.*(((0.5.*(sin(beta).^2).*(cos(theta0).^2) + (cos(beta).^2).*(sin(theta0).^2)).*(cos(gamma).^2)) + (0.5.*(cos(theta0).^2).*(sin(gamma).^2))).*real(at'.*(an-at))) + ((((3/8).*(cos(theta0).^4).*(sin(gamma).^4)) + (((sin(theta0).^4).*(cos(beta).^4) + (3/8).*(cos(theta0).^4).*(sin(beta).^4) + (3/16).*(sin(2.*theta0).^2).*(sin(2.*beta).^2)).*(cos(gamma).^4)) + ((3/16).*((cos(theta0).^4).*(sin(beta).^2) + (sin(2.*theta0).^2).*(cos(beta).^2)).*(sin(2.*gamma).^2))).* ((abs(an-at)).^2)));
f33 = integral2(f3,0,pi/2,0,pi/2);
fvvii2 = C0*C0*n0*(((abs(at)).^2) + f33); % <|fvv(-i,i)|2>
%---------------------------------------------
f4 = @(beta,gamma) (1/(pi/2)).*(1/(pi/2)).*(0.5.*((sin(beta).^2).*(cos(gamma).^2) + (sin(gamma).^2)));
f44 = integral2(f4,0,pi/2,0,pi/2);
fhhii = C0*n0*(at + f44.*(an-at)); % <fhh(+-i,+-i)>
%---------------------------------------------
f5 = @(beta,gamma) (1/(pi/2)).*(1/(pi/2)).*((((sin(beta).^2).*(cos(gamma).^2) + (sin(gamma).^2)).*real(at'.*(an-at))) + ((3/8).*((sin(beta).^4).*(cos(gamma).^4) + 0.5.*(sin(beta).^2).*(sin(2*gamma).^2) + (sin(gamma).^4)).*((abs(an-at)).^2)));
f55 = integral2(f5,0,pi/2,0,pi/2);
fhhii2 = C0*C0*n0*(((abs(at)).^2) + f55); % <|fhh(-i,i)|2>
%---------------------------------------------
f6 = @(beta,gamma) (1/(pi/2)).*(1/(pi/2)).*((((cos(theta0).^2).*(sin(beta).^4) + (sin(2.*beta).^2).*(sin(theta0).^2)).*(cos(gamma).^4)) + ((cos(theta0).^2).*(sin(gamma).^4)) + ((0.5.*(sin(beta).^2).*(cos(theta0).^2) + (cos(beta).^2).*(sin(theta0).^2)).*(sin(2.*gamma).^2)));
f66 = integral2(f6,0,pi/2,0,pi/2);
fvhii2 = C0*C0*(n0/8).*(f66.*((abs(an-at)).^2)); % <|fhh(-i,i)|2>

%---------------------------------------------
sigma0vv = ((2.*pi)./lambda0).*fvvii2.*cos(theta0)./(2.*imag(fvvii));
sigma0hh = ((2.*pi)./lambda0).*fhhii2.*cos(theta0)./(2.*imag(fhhii));
sigma0vh = ((2.*pi)./lambda0).*fvhii2.*cos(theta0)./(imag(fhhii) + imag(fvvii));
sigma0vvdB = 10.*log10(sigma0vv);
sigma0hhdB = 10.*log10(sigma0hh);
sigma0vhdB = 10.*log10(sigma0vh);
end