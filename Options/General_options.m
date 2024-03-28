function pdata = General_options(pdata, varargin)

% The general options includes user-defined values for:
%
%    Water density in [kg/m^3],
%    Gravitational constant in [m/s^2],
%    Water depth in [m],
%    Safety distance for farms in [m],
%    Number of frequencies,
%    Minimum frequency 
%    Maximum frequency,
%    Plot flag,
%    Power conversion chain efficiency,
%    Operational availability,
%    Transmission efficiency,
%    Number of years in the study,
%    Saturation Method, which is set to PowerBound in this release,
%    Wave type (regular or irregular) 
%
% Primary Contributor: Saeed Azad, PhD


% Saturation flag for both F_pto and x
if nargin ==2

    if varargin{1,1}(1,1) == Inf
        pdata.Opt.SatFlag   = 0;
    else
        pdata.Opt.SatFlag   = 1;
    end

    % Type of wave: Default is Regular
    pdata.General.WaveType       = 'Regular';

elseif nargin ==3

    if varargin{1,1}(1,1) == Inf
        pdata.Opt.SatFlag   = 0;
    else
        pdata.Opt.SatFlag   = 1;
    end
    
    % Type of wave: Passed by the user
    pdata.General.WaveType = varargin{1,2};

else

    % Default power saturation is off 
    pdata.Opt.SatFlag   = 0; 

    % Type of wave: Default is Regular 
    pdata.General.WaveType       = 'Regular';

end


% Salt water density [kg/m^3]
pdata.General.rho            = 1025;

% Gravitational constant [m/s^2]
pdata.General.g              = 9.81;

% Depth of the sea [m]
pdata.General.depth          = 50;

% Safety distance [m]
pdata.General.SafeDist       = 10;

% Number of frequencies (Do not change - to be consistent )
pdata.General.nbfreq         = 200;

% Minimum frequency
pdata.General.wMin           = 0.3;

% Maximum frequency
pdata.General.wMax           = 2;

% Flag for plots
pdata.General.plotflag       = 0;

% Power conversion chain efficiency  
pdata.General.PCC_eff        = 0.8;

% Operational availability
pdata.General.Op_av          = 0.95;

% Tranmission efficiency
pdata.General.Trans_eff      = 0.98;


% Number of years of the study
if strcmpi(pdata.General.WaveType,'Irregular')
    pdata.General.n_y        = 30;
elseif strcmpi(pdata.General.WaveType,'Regular')
    pdata.General.n_y        = 1;
end


% Saturation method: only 'PowerBound' for this release
pdata.Opt.SatMethod = 'PowerBound';


if pdata.Opt.SatFlag

    switch pdata.Opt.SatMethod

        case 'PowerBound'

            if nargin == 2

                % Pass user-defined power saturation limit
                pdata.General.PowerLimit   = varargin{1,1}(1,1);

            else

                % Pass default power saturation limit [watt]
                pdata.General.PowerLimit   = 10000;
            end

            % Arbitrary selects a user-defined input to saturate the power
            pdata.General.PowerlimitFlag = 'Arbitrary';

            % Objective and Constraints in the same function (this must be 0 for PowerBound-this release )
            pdata.General.OCFlag         = 0;

    end

end





end
