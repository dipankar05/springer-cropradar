%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 4.7: Relative Dielectric Constant of SOIL

%Description: Code computes the real and imaginary parts of the relative
    %dielectric constant of soil at a given temperature 0<t<40C, frequency,
    %volumetric moisture content, soil bulk density, sand and clay
    %fractions.

%Input Variables:
    %f: frequency in GHz
    %t: temperature in C
    %rho_b: bulk density in g/cm3 (typical value is 1.7 g/cm3)
    %S: Sand Fraction ( 0< S< 1)
    %C: Clay Fraction ( 0< C< 1)
    %mv: Volumetric Water content 0<mv<1

%Output Products:
    %epsr: real part of dielectric constant
    %epsi: imaginary part of dielectric constant
%Book Reference: Section 4-8

%Example call: [A B] = RelDielConst_Soil(f,t,rho_b, mv, S, C)
    %Computes the real and imaginary components of the permitivity of soil
    %based on the temperature value (t) in degrees C and frequency
    %vector (f) and assigns them to vectors A and B respectively

%MATLAB Code

function [epsr epsi] = RelDielConst_Soil(f,t,rho_b, mv,S,C)

f_hz = f * 1.0e9; % transform from GHz to Hz


beta1 = 1.27 - 0.519 * S - 0.152* C;  %eq: 4.68b
beta2 = 2.06 - 0.928 * S - 0.255 * C; %eq: 44.68c 
alpha = 0.65; % eq: 4.68a

eps_0 = 8.854e-12; 

sigma_s = 0;
if f > 1.3
    sigma_s = -1.645 + 1.939 * rho_b - 2.256*S + 1.594 * C; %eq: 4.68d
end
if f >= 0.3 && f <= 1.3 
    sigma_s = 0.0467 + 0.22 * rho_b - 0.411*S + 0.661 *C; %eq: 4.70
end
%Dielectric Constant of Pure Water
    
ew_inf = 4.9; % eq: E.15
ew_0 = 88.045 - 0.4147 * t + 6.295e-4 * t^2 + 1.075e-5 * t^3; %  
tau_w = (1.1109e-10 - 3.824e-12*t +6.938e-14*t^2 - 5.096e-16*t^3)/2/pi; %

epsrW = ew_inf +(ew_0-ew_inf)./(1 + (2*pi*f_hz*tau_w).^2); %

epsiW = 2*pi*tau_w .*f_hz *(ew_0-ew_inf) ./(1 + (2*pi*f_hz*tau_w).^2) + ... 
      (2.65-rho_b)/2.65/mv * sigma_s ./(2*pi*eps_0*f_hz); %

  
% calculating dielectric constant of soil using eq. 4.66a nad 4.66b 
    
    epsr = (1+ 0.66*rho_b + mv^beta1 * epsrW.^alpha - mv).^(1/alpha);
    epsi = mv^beta2 .* epsiW;
    
  end