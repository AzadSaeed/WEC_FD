function [A, B, Fe] = MyMSCoefficients(n_wec,Location_input,wc, Radius, Draft,depth, g, rho, w, varargin)


if length(varargin)>= 1


    % Numerical characteristics
    maxs  = varargin{1,1}.maxs;                % order of interaction
    maxm  = varargin{1,1}.maxm;                % order of series truncation for main fluid
    maxn  = varargin{1,1}.maxn;                % order of series truncation for fluid below body
    maxj  = varargin{1,1}.maxj;                % order of series truncation for diffraction problem
    fmaxw = varargin{1,1}.fmaxw;               % max frequency in rad/s (for wave.m)
    fresw = varargin{1,1}.fresw;               % number of frequencies for wave.m evaluations (linearly spaced from 0 to fmaxw)
    fres  = varargin{1,1}.fres;                % number of frequencies for spectral evaluations (linearly spaced from 0 to fmaxw)

else

    % Numerical characteristics
    maxs  = 5;                % order of interaction
    maxm  = 5;                % order of series truncation for main fluid
    maxn  = 40;                % order of series truncation for fluid below body
    maxj  = 15;                % order of series truncation for diffraction problem
    fmaxw = 2;               % max frequency in rad/s (for wave.m)
    fresw = 200;               % number of frequencies for wave.m evaluations (linearly spaced from 0 to fmaxw)
    fres  = 50;                % number of frequencies for spectral evaluations (linearly spaced from 0 to fmaxw)

end

% Surface, mass, and frequency
S      = pi*(Radius)^2;               % Buoy surface
in.M      = rho*S*Draft;

input_x = Location_input(:,1);
input_y = Location_input(:,2);

wmega_arr = w;
wc        = 1:wc;

% Coefficients to transform to data with dimension
NFA_MS     = rho*pi;                      % Modified normalization factor for A
NFB_MS    = wmega_arr.*rho*pi;           % Modified normalization factor for B
NFFe_MS   = pi*g*rho*(Radius^2);      % Modified normalization factor for Fe


A  = NaN(n_wec,n_wec,length(wc));
B  = NaN(n_wec,n_wec,length(wc));
Fe = NaN(1,n_wec,length(wc));


for i = wc 

        % MS solver: Re-arrange all outputs such that frequency is in the
        % third dimension of the matrix
        wmega            = wmega_arr(i);
        [am,ad,fz]       = wave_s(Radius, Draft, n_wec, input_x', input_y', 0, wmega, depth, g, maxm, maxn, maxj, maxs);
        A(:,:,i)         = am.*NFA_MS;
        B(:,:,i)         = ad*NFB_MS(i);
        Fe(:,:,i)        = fz*NFFe_MS;
        FZ(:,:,i)        = fz;

end

end