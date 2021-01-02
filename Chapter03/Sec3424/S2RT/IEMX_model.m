%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
function sigvh = IEMX_model(fr, sig, L, theta_d, er, sp,xx, auto)

sig = sig * 100 ; % change to cm scale
L = L * 100; % change to cm scale;


%- fr: frequency in GHz
%- sig: rms height of surface in cm
%- L: correlation length of surface in cm
%- theta_d: incidence angle in degrees
%- er: relative permittivity
%- sp: type of surface correlation function

error = 1.0e8;

k = 2*pi *fr/30; % wavenumber in free space. Speed of light is in cm/sec
theta = theta_d .*pi/180; % transform to radian

ks = k * sig; % roughness parameter
kl = k * L;  

ks2 = ks .* ks; 
kl2 = kl.^2;

cs = cos(theta);
s = sin(theta+ 0.001);

s2 = s.^2;

%-- calculation of reflection coefficints
rt = sqrt(er - s2);

rv = (er *cs - rt) ./(er*cs +rt);
rh = (cs - rt)./(cs + rt);

rvh = (rv - rh) ./2;

%-- rms slope values
sig_l = sig/L;
if sp==1                %-- exponential correl func
    rss = sig_l;
end
if sp==2                %-- Gaussian correl func
    rss = sig_l * sqrt(2);
end

if sp==3                 %-- 1.5-power spectra correl func
    rss = sig_l * sqrt(2*xx);
end

%--- Selecting number of spectral components of the surface roughness
if auto == 0 
n_spec = 15; % number of terms to include in the surface roughness spectra
end

if auto == 1
   n_spec = 1; 
   while error > 1.0e-8,
    n_spec = n_spec + 1;
    error = (ks2 .*( 2*cs).^2 ).^n_spec ./ factorial(n_spec); 
   end

end

%-- calculating shadow consideration in single scat (Smith, 1967)

ct = cot(theta+ 0.001);
farg = ct /sqrt(2) ./rss;
gamma = 0.5 *(exp(-farg.^2) / 1.772 / farg - erfc(farg));
Shdw = 1 ./ (1 + gamma);

%-- calculating multiple scattering contribution
%------ a double integration function

svh = dblquad(@(r,phi)xpol_integralfunc(r, phi, sp,xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec), 0.1, 1, 0, pi);

svh = svh *1.e-5; %% un-scale after rescalingin the integrand function.
sigvh = 10*log10(svh .* Shdw);

end

function y = xpol_integralfunc(r, phi, sp, xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec)

cs2 = cs .^2;

r2 = r.^2;
nr = length(r);

sf = sin(phi);
csf = cos(phi);
rx = r .* csf;
ry = r .* sf;

%-- calculation of the field coefficients
rp = 1 + rvh;
rm = 1 - rvh; 

q = sqrt(1.0001 - r2);
qt = sqrt(er - r2);

a = rp ./q;
b = rm ./q;
c = rp ./qt;
d = rm ./qt;

%--calculate cross-pol coefficient
B3 = rx .* ry ./cs;
fvh1 = (b-c).*(1- 3*rvh) - (b - c./er) .* rp; 
fvh2 = (a-d).*(1+ 3*rvh) - (a - d.*er) .* rm;
Fvh = ( abs( (fvh1 + fvh2) .*B3)).^2;


%-- calculate shadowing func for multiple scattering 
au = q ./r ./1.414 ./rss;
fsh = (0.2821./au) .*exp(-au.^2) -0.5 .*(1- erf(au));
sha = 1./(1 + fsh); 

%-- calculate expressions for the surface spectra
wn = spectrm1(sp, xx, kl2, L, rx, ry, s, n_spec, nr);
wm = spectrm2(sp, xx, kl2, L, rx, ry, s, n_spec, nr);


%--compute VH scattering coefficient
acc = exp(-2* ks2 .*cs2) ./(16 .* pi);
vhmnsum = zeros(1,nr);
for n = 1: n_spec
    for m = 1: n_spec
        vhmnsum = vhmnsum + wn(n,:).*wm(m,:) .*(ks2*cs2).^(n+m) ...
            ./factorial(n)./factorial(m); 
    end
end

VH = 4 * acc .* Fvh .* vhmnsum .*r;
y = VH .* sha;

y =y*1.e5; %% rescale so dblquad() works better.
end

function wn = spectrm1(sp, xx, kl2, L, rx, ry, s, np, nr)

wn = zeros(np,nr);

if sp == 1  % exponential
    for n = 1: np       
        wn(n, :) = n* kl2 ./(n.^2 + kl2 *((rx-s).^2+ry.^2)).^1.5;
    end
end

if sp == 2  %  gaussian
    for n = 1: np
        wn(n,:) = 0.5 * kl2 ./n .* exp(-kl2*((rx-s).^2 + ry.^2)/(4*n)) ;
    end
end

if sp== 3 % x-power
    for n = 1: np
        wn(n,:) = kl2 ./(2.^(xx*n-1) .*gamma(xx*n)).* ( ( (rx-s).^2 ...
          + ry.^2).*L).^(xx*n-1) .* besselk(-xx*n+1, L*((rx-s).^2 + ry.^2));
    end
end

end

function wm = spectrm2(sp, xx, kl2, L, rx, ry, s, np, nr)

wm = zeros(np,nr);

if sp == 1  % exponential
    for n = 1: np
        wm(n,:) = n* kl2 ./(n.^2 + kl2 *((rx+s).^2+ry.^2)).^1.5;
    end
end

if sp == 2  %  gaussian
    for n = 1: np
        wm(n,:) = 0.5 * kl2 ./n .* exp(-kl2*((rx+s).^2 + ry.^2)/(4*n)) ;
    end
end

if sp== 3 % x-power
    for n = 1: np
        wm(n,:) = kl2 ./(2.^(xx*n-1) .*gamma(xx*n)).* ( ( (rx+s).^2 ...
          + ry.^2).*L).^(xx*n-1) .* besselk(-xx*n+1, L*((rx+s).^2 + ry.^2));
    end
end
end
