function [sigma0vvdB, sigma0hhdB, sigma0vvdBr, sigma0hhdBr] = Karam_Fung_GRG1988(thetaid)
% Frequency (GHz)
% % High frequency
f = 9.6;
er = 28.04 - 13.34j;

% % Low frequency case
% f = 1.5; % Frequency(GHz);
% er = 28.52 - 2.13j;

%% Vegetation volume backscattering for non-spherical particles
% % Radar incidence angle thi = 30 degree
%thetaid = 30;
thetai = deg2rad(thetaid); 
% In backscatter case, thetas = thetai and phis = phii+pi
dphi = pi;
thetas = thetai;
% In forward scattering case, thetas = pi - thetai and phis = phii

%%
% The scattering amplitudes for a single disc are calculated by using the Generalized Rayleigh-Gans (GRG) approximation
% Consider circular disk like particles
% ------------------------------------------------------------------------------------------------
% Circular disk
% radius = a; thickness = 2t;
a = 5/100; % (5 cm)
t = 0.01/100; % 0.1 mm
% Particle volume
v0 = (4/3)*pi*a*a*t;
m = 2*a/t;
% Calculating demagnetizing factors and modifying functions
% [Ref] Karam, M.A., Fung, A.K. and Antar, Y.M., 1988. Electromagnetic wave scattering from some vegetation samples. IEEE Transactions on Geoscience and Remote Sensing, 26(6), pp.799-808.

gt = (0.5/((m*m)-1))*((((m*m)/sqrt((m*m)-1))*asin(sqrt((m*m)-1)/m))-1);
gn = ((m*m)/((m*m)-1))*(1-((1/sqrt((m*m)-1))*asin(sqrt((m*m)-1)/m)));
at = 1/(((er-1)*gt)+1);
an = 1/(((er-1)*gn)+1);
% modifying function
C = 3*(10e+8); % velocity of light m/sec
lambda = C/(f*(10e+9));
k = (2*pi)/lambda;  % wavenumber
Qsi = k.*sqrt(sin(thetas).^2 + sin(thetai).^2 - (2.*sin(thetas).*sin(thetai).*cos(dphi)));
mus = (2 .* besselj(1,(Qsi.*a)))/(Qsi.*a);

%% Modifying function in case of Rayleigh approiximation
musr = 1;

%%
%% GRG approximation
% Scattering amplitude tensor in backscattered field
Fvvsi  = ((k.*k.*(er-1).*v0)./(4.*pi)) .* ((an.*sin(thetai).*sin(thetas))-(at.*cos(thetai).*cos(thetas).*cos(dphi))) .* mus;
Fhhsi  = ((k.*k.*(er-1).*v0)./(4.*pi)) .* cos(dphi) .* at .* mus;
% Fvhsi  = ((k*k*(er-1)*v0)/(4*pi)) * cos(thetas) * sin(dphi) * at * mus;
% Fhvsi  = ((k*k*(er-1)*v0)/(4*pi)) * cos(thetai) * sin(dphi) * at * mus;

%-----------------------------------------------------------------------------------
sigma0vv = 4.*pi .* (abs(Fvvsi))^2; 
sigma0hh = 4.*pi .* (abs(Fhhsi))^2;
% sigma0hv = 4.*pi .* (abs(Fhvsi))^2;
% sigma0vh = 4.*pi .* (abs(Fvhsi))^2;
%-----------------------------------------------------------------------------------
sigma0vvdB = 10.*log10(sigma0vv);
sigma0hhdB = 10.*log10(sigma0hh);
% sigma0vhdB = 10.*log10(sigma0vh)
% sigma0hvdB = 10.*log10(sigma0hv)

%
%%
%% Rayleigh approximation
% Scattering amplitude tensor in backscattered field
Fvvsi  = ((k.*k.*(er-1).*v0)./(4.*pi)) .* ((an.*sin(thetai).*sin(thetas))-(at.*cos(thetai).*cos(thetas).*cos(dphi))) .* musr;
Fhhsi  = ((k.*k.*(er-1).*v0)./(4.*pi)) .* cos(dphi) .* at .* musr;
% Fvhsi  = ((k*k*(er-1)*v0)/(4*pi)) * cos(thetas) * sin(dphi) * at * mus;
% Fhvsi  = ((k*k*(er-1)*v0)/(4*pi)) * cos(thetai) * sin(dphi) * at * mus;

%-----------------------------------------------------------------------------------
sigma0vv = 4.*pi .* (abs(Fvvsi))^2; 
sigma0hh = 4.*pi .* (abs(Fhhsi))^2;
% sigma0hv = 4.*pi .* (abs(Fhvsi))^2;
% sigma0vh = 4.*pi .* (abs(Fvhsi))^2;
%-----------------------------------------------------------------------------------
sigma0vvdBr = 10.*log10(sigma0vv);
sigma0hhdBr = 10.*log10(sigma0hh);
% sigma0vhdB = 10.*log10(sigma0vh)
% sigma0hvdB = 10.*log10(sigma0hv)

end