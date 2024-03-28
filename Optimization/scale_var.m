function [y,D_scale,b_scale] = scale_var(x,l,u)

% Scaling function to obtain range and offset
% Inputs
% x ------- main variable
% l ------- Lower bound of the variable
% u ------- Upper bound of the variable
%
% Oututs
% y       ------- scaled variabl
% D_scale ------- Diagonal matrix
% b_scale ------- Offset 

    x = x(:); l = l(:); u = u(:);

    d = 2./(u-l);
    D_scale = diag(d);
    b_scale = -(u+l)./(u-l);
    y = D_scale*x + b_scale;

end