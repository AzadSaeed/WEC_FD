function pdata = PC_Opt_options(pdata)

% Surrogateopt, fmincon, GA
if ~isfield(pdata.Opt,'Solver')
    pdata.Opt.Solver      = 'GA';
end

% Max function evaluation
pdata.Opt.MaxFuncEval    = 500;


% Number of plant variables [Radius (m), Radius/Draft]
pdata.Opt.np             = 2;

% Number of location variables [x (m), y (m)], first array at [0,0]
pdata.Opt.nl             = 0;



% Define bounds on plant variables
pdata.Opt.lbp                = [pdata.Opt.R_min; pdata.Opt.RD_min];
pdata.Opt.ubp                = [pdata.Opt.R_max; pdata.Opt.RD_max];



% Number of Control Variables [Damping (N.s/m), stiffness (N/m)]
[lbc,ubc, nc]   = Control_Bounds(pdata.General.n_wec,pdata.Opt.B_PTO_min,...
    pdata.Opt.k_PTO_min,pdata.Opt.B_PTO_max,pdata.Opt.k_PTO_max,...
    pdata.Opt.CtrlFlag);
pdata.Opt.lbc = lbc;
pdata.Opt.ubc = ubc;
pdata.Opt.nc  = nc;

% Number of total variables
pdata.Opt.nt                 = pdata.Opt.nc+pdata.Opt.np+pdata.Opt.nl;

% Assemble all bounds
pdata.Opt.lb                 = [pdata.Opt.lbc; pdata.Opt.lbp];
pdata.Opt.ub                 = [pdata.Opt.ubc; pdata.Opt.ubp];

% Prescribe WEC Locations
if ~isfield(pdata.General,'Location_Flag')
    Location_Flag = 'FAR-SYMMETRICAL';
    pdata.General.Location_Flag = Location_Flag;
end


% Prescribe Location
[x,y]       = Locations(pdata.General.Location_Flag,pdata.General.n_wec);
pdata.Opt.x = x;
pdata.Opt.y = y;

end