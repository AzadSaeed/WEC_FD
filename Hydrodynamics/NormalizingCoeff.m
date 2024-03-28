function [A_cf, B_cf, Fe_cf] = NormalizingCoeff(Radius, Draft, Omega, rho, g) 

Omega     = Omega(:);

A_cf      = rho*pi*(Radius.^3);                      % New modified normalization factor for A
B_cf      = Omega.*rho*pi.*(Radius.^3);              % New modified normalization factor for B
Fe_cf     = (rho*pi.*(Radius.^2).*(Draft))*g;         % Conventional normalization factor for Fe

end