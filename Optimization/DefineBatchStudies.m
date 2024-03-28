function DefineBatchStudies(pdata)
% This function allows the user to define various cases to study in batch.
% All the requires definitions are created in cell arrays, and are passed to
% the RunBatchStudies function to be used in for-loops for batch implementation.
% 
%
% Saeed Azad, PhD


CaseStudy1 = 'P_opt';    % Plant-only optimization
CaseStudy2 = 'C_opt';    % Control-only optimization
CaseStudy3 = 'L_opt';    % Layout-only optimization
CaseStudy4 = 'LP_opt';   % Concurrent Layout and plant optimization
CaseStudy5 = 'PC_opt';   % Concurrent plant and control optimization
CaseStudy6 = 'PCL_opt';  % Concurrent plant, control, and layout optimization
CaseStudy7 = 'Sim';      % Simulations

% Run an individual study

% Case study:'P_opt','C_opt','L_opt','LP_opt','PC_opt','PCL_opt', 'Sim'
Studies    = {CaseStudy1,CaseStudy2};

% Regions = 'AlaskaCoasts','EastCoast','PacificIslands','WestCoast'
Regions    = {'WestCoast','EastCoast'};

% Hedrodynamics:
%  - Multi scattering 'MS'
%  - Surrogate models with Many-body exapnsion 'SM_MBE'
%  - Multi scattering with Many-body expansion 'MS-MBE'
Hydro      = {'MS','SM_MBE'};

% Controller flag:
% Device level: 'Individual'
% farm level: 'farm'
Ctrl       = {'Farm'};

% Starting rng for optimization studies
IP         = {15,10};

% Prescribed layouts (must be defined)
Location_Flag = {'NEAR-SYMMETRICAL', 'Far-symmetrical'};

% Optimization solver
Solver     = {'fmincon'};

% Power saturation limit in Watts, use INF for no imit
PowerLimit    = {1000, 2000};

% Wave type
WaveType = {'Regular', 'Irregular'};

% Run Batch Studies
pdata = RunBatchStudies(Regions,Studies,Hydro,IP,Ctrl,Location_Flag,...
    Solver,PowerLimit, WaveType, pdata);

end