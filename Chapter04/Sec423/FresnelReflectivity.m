function [gammav, gammah] = Fresn_Refl(eps,theta)
% calculates Fresnel reflectivities of v and h-polarizations at given set
% of incidence angles.

[rho_v, rho_h] = refl_coef(theta, 1, eps);
gammav = (abs(rho_v)).^2;
gammah = (abs(rho_h)).^2;

end

%----------------------------------------------------------------
function [rho_v, rho_h] = refl_coef(the1, eps1, eps2)

% calculates the v and h-polarized reflection coefficients of a plane
% dielectric surface

n1 = sqrt(eps1);
n2 = sqrt(eps2);
costh2 = sqrt(1 - (n1.*sin(the1)./n2).^2);


rho_v = -(n2.*cos(the1) - n1.*costh2)./(n2.*cos(the1)+n1.*costh2);
rho_h = (n1.*cos(the1) - n2.*costh2)./(n1.*cos(the1)+n2.*costh2);
end
