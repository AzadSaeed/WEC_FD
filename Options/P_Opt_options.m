function pdata = P_Opt_options(pdata)

% Surrogateopt, fmincon, GA
if ~isfield(pdata.Opt,'Solver')
    pdata.Opt.Solver      = 'fmincon';
end


% Number of plant variables [Radius (m), Radius/Draft]
pdata.Opt.np             = 2;

% Number of location variables [x (m), y (m)], first array at [0,0]
pdata.Opt.nl             = 0;

% Number of Control Variables [Damping (N.s/m), stiffness (N/m)]
pdata.Opt.nc             = 0;

% Number of total variables
pdata.Opt.nt             = pdata.Opt.nc+pdata.Opt.np+pdata.Opt.nl;

% Define bounds on plant variables
pdata.Opt.lbp            = [pdata.Opt.R_min; pdata.Opt.RD_min];
pdata.Opt.ubp            = [pdata.Opt.R_max; pdata.Opt.RD_max];

% Assemble all bounds
pdata.Opt.lb             = [pdata.Opt.lbp];
pdata.Opt.ub             = [pdata.Opt.ubp];

% Prescribe WEC Locations
if ~isfield(pdata.General,'Location_Flag')
    Location_Flag = 'FAR-SYMMETRICAL';
    pdata.General.Location_Flag = Location_Flag;
end

% Prescribe Location
[x,y]                    = Locations(pdata.General.Location_Flag,pdata.General.n_wec);
pdata.Opt.x              = x;
pdata.Opt.y              = y;

% PTO stiffness
pdata.Opt.k_pto          = -5*10^3*ones(1,pdata.General.n_wec);

% PTO damping
pdata.Opt.B_pto          = pdata.Opt.B_PTO_max*ones(1,pdata.General.n_wec);

end