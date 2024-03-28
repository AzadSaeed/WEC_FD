function pdata = C_Opt_options(pdata)

% Surrogateopt, fmincon, GA
if ~isfield(pdata.Opt,'Solver')
    pdata.Opt.Solver      = 'fmincon';
end

pdata.Opt.MaxFuncEval    = 500;

% Number of plant optimization variables
pdata.Opt.np             = 0;

% Number of layout optimization variables
pdata.Opt.nl             = 0;


try
    % Load Radius and slenderness ratio from a p_opt solution
    sol_p = load('Solution_xx');
    pdata.Opt.R          = sol_p.pdata.Results.Radius;
    pdata.Opt.RD         = sol_p.pdata.Results.RD;
    clear sol_p

catch

    % default radius and slenderness ratios
    pdata.Opt.R          = 5;
    pdata.Opt.RD         = 5;

end


if ~isfield(pdata.General,'Location_Flag')
    Location_Flag = 'FAR-SYMMETRICAL';
    pdata.General.Location_Flag = Location_Flag;
end

% Prescribe Location
[x,y]       = Locations(pdata.General.Location_Flag,pdata.General.n_wec);
pdata.Opt.x = x;
pdata.Opt.y = y;


% Number of Control Variables [Damping (N.s/m), stiffness (N/m)]
[lbc,ubc, nc]   = Control_Bounds(pdata.General.n_wec,pdata.Opt.B_PTO_min,...
    pdata.Opt.k_PTO_min,pdata.Opt.B_PTO_max,pdata.Opt.k_PTO_max,...
    pdata.Opt.CtrlFlag);
pdata.Opt.lbc = lbc;
pdata.Opt.ubc = ubc;
pdata.Opt.nc = nc;


switch upper(pdata.Opt.CtrlFlag)

    case 'INDIVIDUAL'

        % Optimize control of each individual device
        pdata.Opt.nc   = 2*pdata.General.n_wec;
        pdata.Opt.lbc  = repmat([pdata.Opt.B_PTO_min; pdata.Opt.k_PTO_min],[pdata.General.n_wec,1]);
        pdata.Opt.ubc  = repmat([pdata.Opt.B_PTO_max; pdata.Opt.k_PTO_max],[pdata.General.n_wec,1]);

    case 'FARM'

        % Optimize one control for the entire farm
        pdata.Opt.nc   = 2;
        pdata.Opt.lbc  = [pdata.Opt.B_PTO_min ; pdata.Opt.k_PTO_min];
        pdata.Opt.ubc  = [pdata.Opt.B_PTO_max ; pdata.Opt.k_PTO_max];
end


% Number of total variables
pdata.Opt.nt                 = pdata.Opt.nc+pdata.Opt.np+pdata.Opt.nl;

% Assemble all bounds
pdata.Opt.lb                 = [pdata.Opt.lbc];
pdata.Opt.ub                 = [pdata.Opt.ubc];

end