function pdata = PCL_Opt_options(pdata, varargin)

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
pdata.Opt.MaxFuncEval    = 50;

% Number of plant variables [Radius (m), Radius/Draft]
pdata.Opt.np             = 2;

% Number of location variables [x (m), y (m)], first array at [0,0]
pdata.Opt.nl             = 2*(pdata.General.n_wec-1);


% Number of Control Variables [Damping (N.s/m), stiffness (N/m)]
[lbc,ubc, nc]   = Control_Bounds(pdata.General.n_wec,pdata.Opt.B_PTO_min,...
    pdata.Opt.k_PTO_min,pdata.Opt.B_PTO_max,pdata.Opt.k_PTO_max,...
    pdata.Opt.CtrlFlag);
pdata.Opt.lbc = lbc;
pdata.Opt.ubc = ubc;
pdata.Opt.nc  = nc;


if pdata.Opt.Hybrid

    % Load from a previous solution
    filename = strcat('Solution_GA_PCL_opt_50_SM_MBE_Farm_rng_15_Plim_None');

    sol_A = load(strcat(pdata.General.pathstr,filesep,'Solution',filesep,...
        pdata.General.CaseStudy, filesep,'5_WEC',filesep,'Irregular',filesep, 'Results',...
       filesep, Region, filesep,filename));

    pdata.Opt.R     = sol_A.pdata.Results.Radius;
    pdata.Opt.RD    = sol_A.pdata.Results.RD;

    plant_ = [pdata.Opt.R;pdata.Opt.RD];

    pdata.Opt.k_pto = sol_A.pdata.Results.PTO_stiffness';
    pdata.Opt.B_pto = sol_A.pdata.Results.PTO_damping';

    if strcmpi(pdata.Opt.CtrlFlag,'FARM')
        control_ = [pdata.Opt.B_pto(1,1); pdata.Opt.k_pto(1,1)];
    elseif strcmpi(pdata.Opt.CtrlFlag,'INDIVIDUAL')
        for ii = 1:pdata.General.n_wec
            control_(ii,:) = [pdata.Opt.B_pto(ii,1);pdata.Opt.k_pto(ii,1)];
        end
    end

    pdata.Opt.x0    = sol_A.pdata.Results.x;
    pdata.Opt.y0    = sol_A.pdata.Results.y;
    loc_ = [];
    for ii = 1:pdata.General.n_wec
        loc_((2*ii)-1:2*(ii),:) = [pdata.Opt.x0(ii,1);pdata.Opt.y0(ii,1)];
    end
    loc_ = loc_(3:end,:);

    clear sol_A

    % Define the Stating point
    pdata.Opt.X_Init = [control_;plant_;loc_];

    % Always use fmincon with multi scattering for hybrid optimization
    pdata.Opt.Solver      = 'fmincon';
    pdata.HydroFlag          = 'MS';

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
pdata.Opt.lb                 = [pdata.Opt.lbc; pdata.Opt.lbp; pdata.Opt.lbl];
pdata.Opt.ub                 = [pdata.Opt.ubc; pdata.Opt.ubp; pdata.Opt.ubl];



end