% This is the main file that allows the user to run different simulations
% or optimization problems within the context of the tool.
%
%
% While some problem options, such as the desired case study, and the
% number of WEC devices to study can be set on this file, the majority of
% the settings and problems options can be specified and updated by the
% user through ProblemOptions.m
%
% The toolbox enable the optimization of WEC plant (radius, slenderness ratio)
% control (PTO damping, PTO stiffness), and layout (x_wecs, y_wecs).
% Individual and integrated optimization studies are enables in this
% toolbox. Therefore, the user can run any of the problem combinations possible,
% including plant optimization ('P_opt'), control optimizatio ('C_opt')  
% 
% This work, its creation, and results were made possible through the
% financial support from National Science Foundation Engineering, Design
% and Systems Engineering Program, USA under grant number CMMI-2034040
%
%
% Primary contributor: Saeed Azad, PhD


clear; clc; close all;


% Establish the current directory
pdata = EstablishDirectory;


% List of all Case Studies:
CaseStudy1 = 'P_opt';    % Plant-only optimization
CaseStudy2 = 'C_opt';    % Control-only optimization
CaseStudy3 = 'L_opt';    % Layout-only optimization
CaseStudy4 = 'LP_opt';   % Concurrent Layout and plant optimization
CaseStudy5 = 'PC_opt';   % Concurrent plant and control optimization
CaseStudy6 = 'PCL_opt';  % Concurrent plant, control, and layout optimization
CaseStudy7 = 'Sim';      % Simulations


% Select the number of WECs in the farm (minimum is 1)
pdata.General.n_wec      = 5;


% Select if you want to run a single study (0) or a batch of studies (1)
BatchStudies = 0;


% Run an individual study
if BatchStudies== 0

    % Case study:'P_opt','C_opt','L_opt','LP_opt','PC_opt','PCL_opt', 'Sim'
    Studies    = CaseStudy1;
    pdata.General.CaseStudy  = Studies;

    if strcmpi(pdata.General.CaseStudy,'Sim')


        % Define simulation conditions
        pdata = Sim_Func(pdata);

        % Process simulation cases and results
        pdata = Process_Sim(pdata) ;


    else % All other optimization cases

        % Regions = 'AlaskaCoasts','EastCoast','PacificIslands','WestCoast'
        Regions    = 'WestCoast';


        % Hedrodynamics:
        %  - Multi scattering 'MS'
        %  - Surrogate models with Many-body exapnsion 'SM_MBE'
        %  - Multi scattering with Many-body expansion 'MS-MBE'
        Hydro      = 'SM_MBE';

        % Controller flag:
        % Device level: 'Individual'
        % farm level: 'farm'
        Ctrl       = 'Farm';

        % Starting rng for optimization studies
        IP         = 15;

        % Prescribed layouts (must be defined)
        Location_Flag = {'NEAR-SYMMETRICAL'};

        % Optimization solver
        Solver     = {'fmincon'};

        % Power saturation limit in Watts, use INF for no imit
        PowerLimit    = Inf;


        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,Regions,Hydro,Ctrl,IP, Location_Flag{1,1},Solver{1,1},PowerLimit, 'Regular');

        % Load all surrogate models once-if needed
        if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
            % 'CMT' or 'Non-CMT'
            pdata = LoadALLSMs(pdata,'CMT');
        end

        pdata = RunCases(pdata);

    end


elseif BatchStudies == 1

    % Define new batch Studies
    DefineBatchStudies(pdata)

    % Or run ACC cases
    % ACC_Runs(pdata)
end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function pdata = EstablishDirectory

% This function establishes the current directory of the project, in order
% to later create solution direcotries in the right location
%
%
% Contributor: Saeed Azad, PhD

% Full directory name
pdata.General.dir_path = mfilename('fullpath');

% Directory path
pdata.General.pathstr  = fileparts(pdata.General.dir_path);

% Go to the current directory
cd(pdata.General.pathstr)

% Add the entire directory to MATLAB's path
addpath(genpath(pdata.General.pathstr))

end





