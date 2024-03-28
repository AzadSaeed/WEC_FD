function pdata = Optimization_options(pdata,varargin)

% This function defines various optimization-related parameters, options,
% and settings. In order, the varargin arguments can take the following
% parameters: 
% 
% Region, Hydrodynamic Flag, Controller Flag, Initial Optimization RNG 
% value, Prescribed layout, Optimization solver, and Power Saturation Limit
%
% 
% The options defined in this function include:
%
%   HydroFlag ----- Method for calculating hydrodynamic coefficients.
%                   current options are 'MS' for multi scattering, 'SM_MBE' 
%                   for surrogate modeling with many-body expansion, and
%                   'MS_MBE' multiscattring with many-body expansion
%
%   CtrlFlag ----- Flag for controller: farm level or device level. 'Farm'
%                  level specifies control parameters uniformly for all of
%                  the devices across the farm, while 'individual' entails
%                  the determination of control parameters for each device 
%
%   rngval ------- RNG value for random yet repeatable initialization of
%                  gradient-based optimization algorithm
% 
%   Loc_Flag ----- Name of the layout prescribed to run the case study. The
%                  fact that layout is prescribed means that this flag is
%                  only consequential for cases where layout is not
%                  optimized. The available prescribed layouts can be found
%                  in the ../options/Locations. Any new layout can be added
%                  to the Locations.m function
%
%   Solver ------- Solver is a flag to let the toolbox know which
%                  optimization solver should be used. Currently the
%                  toolbox can take 'fmincon', 'GA', and 'surrogateopt'
%
%   ObjForm ------ Flag for objective function definition. Cane take 'power'
%                  or 'powerpervolume'
%
%   Parflag  ----- Flag for parallel computaiton. You can assign the number
%                  of cores to be used by changing NumWorkers 
%
%  Hybrid -------- Flag for hybrid optimization
%
%  Bounds -------- Includes minimum and maximum values for plant, control,
%                  and layout and draft constraint. The plant optimization  
%
%  Case-specific options are defined in their associated functions as
%  
%  Solver options are also defined in function. These include maximum
%  number of function evaluaitons in fmincon, population number in GA, etc.
%
% Primary Contributor: Saeed Azad, PhD

if nargin ==1

    % Default values
    
    % HydrodynamicFlag 
    HydroFlag  = 'MS';

    % Controller Flag
    CtrlFlag   = 'FARM';

    % Initial Optimization RNG value
    rngval     = 15;

    % prescribed layout name
    Loc_Flag   = 'FAR-SYMMETRICAL';

    % Optimization solver
    Solver     = 'fmincon';

elseif nargin == 2

    if length(varargin{1,1}) == 1

        Region = varargin{1,1}{1,1};
        HydroFlag = 'MS';
        CtrlFlag = 'FARM';
        rngval   = 15;
        Loc_Flag = 'FAR-SYMMETRICAL';
        Solver    = 'fmincon';


    elseif length(varargin{1,1}) == 2

        Region = varargin{1,1}{1,1};
        HydroFlag = varargin{1,1}{1,2};
        CtrlFlag = 'FARM';
        rngval   = 15;
        Loc_Flag = 'FAR-SYMMETRICAL';
        Solver    = 'fmincon';

    elseif length(varargin{1,1}) == 3

        Region = varargin{1,1}{1,1};
        HydroFlag = varargin{1,1}{1,2};
        CtrlFlag = varargin{1,1}{1,3};
        rngval   = 15;
        Loc_Flag = 'FAR-SYMMETRICAL';
        Solver    = 'fmincon';

    elseif length(varargin{1,1}) == 4
        
        Region = varargin{1,1}{1,1};
        HydroFlag = varargin{1,1}{1,2};
        CtrlFlag = varargin{1,1}{1,3};
        rngval   = varargin{1,1}{1,4};
        Loc_Flag = 'FAR-SYMMETRICAL';
        Solver    = 'fmincon';

   elseif length(varargin{1,1}) == 5

        Region = varargin{1,1}{1,1};
        HydroFlag = varargin{1,1}{1,2};
        CtrlFlag = varargin{1,1}{1,3};
        rngval   = varargin{1,1}{1,4};
        Loc_Flag = varargin{1,1}{1,5};
        Solver    = 'fmincon';

    elseif length(varargin{1,1}) == 6

        Region = varargin{1,1}{1,1};
        HydroFlag = varargin{1,1}{1,2};
        CtrlFlag = varargin{1,1}{1,3};
        rngval   = varargin{1,1}{1,4};
        Loc_Flag = varargin{1,1}{1,5};
        Solver    = varargin{1,1}{1,6};

    elseif length(varargin{1,1}) == 7

        Region = varargin{1,1}{1,1};
        HydroFlag = varargin{1,1}{1,2};
        CtrlFlag = varargin{1,1}{1,3};
        rngval   = varargin{1,1}{1,4};
        Loc_Flag = varargin{1,1}{1,5};
        Solver    = varargin{1,1}{1,6};
        PowerLimit = varargin{1,1}{1,7};
        pdata.General.PowerLimit = PowerLimit;

     elseif length(varargin{1,1}) == 8

        Region = varargin{1,1}{1,1};
        HydroFlag = varargin{1,1}{1,2};
        CtrlFlag = varargin{1,1}{1,3};
        rngval   = varargin{1,1}{1,4};
        Loc_Flag = varargin{1,1}{1,5};
        Solver    = varargin{1,1}{1,6};
        PowerLimit = varargin{1,1}{1,7};
        pdata.General.PowerLimit = PowerLimit;

    end
end

% Pass argument
pdata.General.Location_Flag = Loc_Flag;
pdata.Opt.Solver         = Solver;

% Objective funciton power, or PowerPerVolume
pdata.Opt.ObjForm        = 'PowerPerVolume';

% Parallel flag
pdata.Opt.Parflag        = false;

% If parallel, how many cores? 
if pdata.Opt.Parflag
    feature('numcores');
    c = parcluster('local');
    c.NumWorkers;
    c.NumWorkers = 32;
    saveProfile(c);
    parpool('local',c.NumWorkers);
end


% Hydrodynamics flag: 
% Multi scattring 'MS' 
% Surrogate Modeling suing Many-Body Expansion 'SM_MBE'
% Multi scattering with Many-Body Expansion 'MS_MBE'
pdata.HydroFlag          = HydroFlag;


% Control flag: individual or farm
pdata.Opt.CtrlFlag       = CtrlFlag;


% Area flag: half or full: half only considers the rigth plane, while full
% considers both righ and left planes
pdata.Opt.AreaFlag       = 'half';

% Hybrid optimization flag: Turn this option on if you want to use a 
% gradient-based optimizer such as fmincon to start with a solution 
% obtained by genetic algorithm. Similarly this can be used to start with a
% solution obtained using SM_MBE and then refined using MS 
pdata.Opt.Hybrid         = 0;


% Rng value for the starting point 
pdata.Opt.rng            = rngval;


% CheckPoint Flag for saving time -  computationally expensive, so keep 0
pdata.Opt.CheckPointFlag = 0;

% Radius [m]
pdata.Opt.R_min          = 0.5;                   
pdata.Opt.R_max          = 10;                      

% Radius to draft ratio  [m] 
pdata.Opt.RD_min         = 0.2;                   
pdata.Opt.RD_max         = 10;                    

% Draft [m] -- defined to be half of the overall height
pdata.Opt.D_min          = 0.5;                   
pdata.Opt.D_max          = 20;                

% Distance [m]
% pdata.Opt.Dist_min       = 2*pdata.Opt.R_min+pdata.General.SafeDist;      
% pdata.Opt.Dist_max       = 1000;

% Angle
% pdata.Opt.Angle_min      = 0;
% pdata.Opt.Angle_max      = pi;

% PTO [Damping (N.s/m), stiffness (N/m)]
pdata.Opt.k_PTO_min      = -5*10^5;
pdata.Opt.k_PTO_max      = 5*10^5;
pdata.Opt.B_PTO_min      = 0;
pdata.Opt.B_PTO_max      = 5*10^5;

% For regular waves, turn FindOptimalPTO on if you would like to find to
% analytically find the optimal PTO solution.
if strcmpi(pdata.General.WaveType, 'regular')
    pdata.General.FindOptimalPTO = 0;
end


% Case-specific options
switch upper(pdata.General.CaseStudy)


    case 'P_OPT'

        pdata = P_Opt_options(pdata);

    case 'C_OPT'

        pdata = C_Opt_options(pdata);

    case 'L_OPT'

        pdata = L_Opt_options(pdata,Region);

    case 'LP_OPT'

        pdata = LP_Opt_options(pdata,Region);

    case 'PC_OPT'

        pdata = PC_Opt_options(pdata);

    case "PCL_OPT"

        pdata = PCL_Opt_options(pdata,Region);
        

end


% Solver options
switch upper(pdata.Opt.Solver)

    case 'FMINCON'

        pdata.Opt.MaxFuncEval    = 500;
        pdata.Opt.Tolerance      = 10^-8;
        pdata.Opt.FiniteDifferenceStepSize = 10^-6;
        pdata.Opt.Algorithm      = 'interior-point';
        pdata.Opt.Display        = 'Iter';

    case 'GA'

        % Population Size (default'50 when numberOfVariables <= 5, else 200')
        if pdata.Opt.nt <= 5
            pdata.Opt.PS = 50;
        else
            pdata.Opt.PS = 200;
        end

        % MaxGeneration (default = 100*pdata.Opt.nt)
        pdata.Opt.MG = 50*pdata.Opt.nt;
        pdata.Opt.ConstraintTolerance = 1e-6;
        pdata.Opt.CrossoverFraction = 0.8;
        pdata.Opt.FunctionTolerance = 10^-16;
        pdata.Opt.MaxStallGen = 50;
        
        % Define selection function (tournament) with size 2
        pdata.Opt.Se_Func = {@selectiontournament,2};

    case 'SURROGATEOPT'

        pdata.Opt.PlotFunc = 'surrogateoptplot';


end




end