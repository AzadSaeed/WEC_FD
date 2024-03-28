function pdata = L_Opt_options(pdata,varargin)


if nargin == 1

    Region = '';

elseif nargin == 2
    
    Region = varargin{1,1};

end


% Surrogateopt, fmincon, GA
if ~isfield(pdata.Opt,'Solver')
    pdata.Opt.Solver      = 'GA';
end


pdata.Opt.MaxFuncEval    = 1500;

% Number of plant optimization variables
pdata.Opt.np             = 0;

% Number of layout optimization variables
pdata.Opt.nl             = 2*(pdata.General.n_wec-1);

% Number of Control Variables [Damping (N.s/m), stiffness (N/m)]
pdata.Opt.nc             = 0;


if pdata.Opt.Hybrid

    % Load from a previous solution 
    filename = strcat('Solution_surrogateopt_L_opt_1_MS_rng_15_Plim_None');

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
    pdata.Opt.x0    = sol_A.pdata.Results.x;
    pdata.Opt.y0    = sol_A.pdata.Results.y;
    loc_ = [];
    for ii = 1:pdata.General.n_wec
        loc_((2*ii)-1:2*(ii),:) = [pdata.Opt.x0(ii,1);pdata.Opt.y0(ii,1)];
    end
    loc_ = loc_(3:end,:);

    clear sol_A

    % Define the Stating point
    pdata.Opt.X_Init = [loc_];

    % Always use fmincon with multi scattering for hybrid optimization
    pdata.Opt.Solver      = 'fmincon';
    pdata.HydroFlag          = 'MS';


else

    % Radius of WEC
    pdata.Opt.R  = 0.5;

    % WEC slenderness ratio 
    pdata.Opt.RD = 1;

    % PTO stiffness 
    pdata.Opt.k_pto = -5*10^3*ones(1,pdata.General.n_wec);

    % PTO damping 
    pdata.Opt.B_pto = 5*10^5*ones(1,pdata.General.n_wec);

end


% Number of total variables
pdata.Opt.nt   = pdata.Opt.nc+pdata.Opt.np+pdata.Opt.nl;


% Define bounds on farm area
[lbl, ubl]     = Area_Bounds(pdata.General.n_wec,pdata.Opt.AreaFlag);
pdata.Opt.lbl  = lbl;
pdata.Opt.ubl  = ubl;


% Assemble all bounds
pdata.Opt.lb   = [pdata.Opt.lbl];
pdata.Opt.ub   = [pdata.Opt.ubl];

end