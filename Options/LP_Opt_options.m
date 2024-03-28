function pdata = LP_Opt_options(pdata, varargin)

if nargin == 1

    Region = '';

elseif nargin == 2
    
    Region = varargin{1,1};

end

% Surrogateopt, fmincon, GA
if ~isfield(pdata.Opt,'Solver')
    pdata.Opt.Solver      = 'GA';
end

% Max function evaluation
pdata.Opt.MaxFuncEval    = 1500;

% Number of plant variables [Radius (m), Radius/Draft]
pdata.Opt.np             = 2;

% Number of location variables [x (m), y (m)], first array at [0,0]
pdata.Opt.nl             = 2*(pdata.General.n_wec-1);


% Number of Control Variables [Damping (N.s/m), stiffness (N/m)]
pdata.Opt.nc             = 0;


if pdata.Opt.Hybrid

    % Load from a previous solution 
    filename = strcat('Solution_GA_LP_opt_1_SM_MBE_rng_15_Plim_None');

    sol_A = load(strcat(pdata.General.pathstr,filesep,'Solution',filesep,...
        pdata.General.CaseStudy, filesep,'5_WEC',filesep,'Irregular',filesep, 'Results', filesep,...
        Region, filesep,filename));

    pdata.Opt.R     = sol_A.pdata.Results.Radius;
    pdata.Opt.RD    = sol_A.pdata.Results.RD;

    pdata.Opt.k_pto = sol_A.pdata.Results.PTO_stiffness;
    pdata.Opt.B_pto = sol_A.pdata.Results.PTO_damping;

    pdata.Opt.x0    = sol_A.pdata.Results.x;
    pdata.Opt.y0    = sol_A.pdata.Results.y;

    % Assemble initial point
    plant_ = [pdata.Opt.R;pdata.Opt.RD];

    pdata.Opt.x0    = sol_A.pdata.Results.x;
    pdata.Opt.y0    = sol_A.pdata.Results.y;
    loc_ = [];
    for ii = 1:pdata.General.n_wec
        loc_((2*ii)-1:2*(ii),:) = [pdata.Opt.x0(ii,1);pdata.Opt.y0(ii,1)];
    end
    loc_ = loc_(3:end,:);

    clear sol_A

    % Define the Stating point
    pdata.Opt.X_Init = [plant_;loc_];

    % Always use fmincon with multi scattering for hybrid optimization
    pdata.Opt.Solver      = 'fmincon';
    pdata.HydroFlag          = 'MS';

else

    % PTO stiffness - change to optimal plant from 'C_OPT' case
    pdata.Opt.k_pto          = -5*(10^3)*ones(1,pdata.General.n_wec);

    % PTO damping - change to optimal plant from 'C_OPT' case
    pdata.Opt.B_pto          = 5*(10^5)*ones(1,pdata.General.n_wec);

end



% Number of total variables
pdata.Opt.nt                 = pdata.Opt.nc+pdata.Opt.np+pdata.Opt.nl;


% Define bounds on plant variables
pdata.Opt.lbp                = [pdata.Opt.R_min; pdata.Opt.RD_min];
pdata.Opt.ubp                = [pdata.Opt.R_max; pdata.Opt.RD_max];

% Define bounds on farm area
[lbl, ubl]                   = Area_Bounds(pdata.General.n_wec,pdata.Opt.AreaFlag);
pdata.Opt.lbl                = lbl;
pdata.Opt.ubl                = ubl;

% Assemble all bounds
pdata.Opt.lb                 = [pdata.Opt.lbp; pdata.Opt.lbl];
pdata.Opt.ub                 = [pdata.Opt.ubp; pdata.Opt.ubl];

end