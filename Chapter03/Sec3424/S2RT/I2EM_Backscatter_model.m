%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 10.1: I2EM Backscattering from Single-Scale Random Surface

%Description: Code computes sigma_0_vv, sigma_0_hh, and
%sigma_0_hv for single-scale random surface with a specified correlation
%function

%Input Variables: 
    %er: complex dielectric constant of the scattering medium
    %thi: Incidence angle (deg)
    
    %sp: type of correlation function: 1- exponential, 2- Gaussian, 3-
            %x-power
    %xx: coefficient for the x-power correlation function 
    %sig: rms height (m)
    %L: correlation length (m)
    %fr: frequency (GHz)
    
%Output Products:
    % sigma_0_vv, sigma_0_hh, and sigma_0_hv in dB

%Book Reference: Section 10-3.9

%Matlab Code: 

function [sigma_0_vv sigma_0_hh sigma_0_hv] = I2EM_Backscatter_model(...
    fr, sig, L, thi, er, sp, xx)


%--the code calls two other functions. The first is the I2EM_Bistatic_model
%code which is used to calculate the co-polarized responses in backscatter.
%The second is the IEMX_model used to calculate the cross-polarized
%response. The later code is slow as it involved double integrations.


%--- The co-pol components:
ths = thi;
phs = 180;
[sigma_0_vv sigma_0_hh] = I2EM_Bistat_model(fr, sig, L, thi, ths, phs, er, sp, xx);


%--- The cross-pol component
auto = 1; % auto = 1 allows for automatic selection of the number of spectral components 
          % auto = 0 forces the number of spectral components to 15 always.
          % selection of auto = 1 results in a slower code

[sigma_0_hv] = IEMX_model(fr, sig, L, thi, er, sp,xx, auto);
 
end






