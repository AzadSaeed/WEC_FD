% This script is written to test the hybrid option in the code. The hybrid
% option uses a solution obtained by a global optimization solver such as
% GA, using surrogate models, as the starting point for a gradient-based 
% solver using multi-scattring.
% Since this function is mainly useful for studies involving layout
% optimization, it is only available for 'L_opt', 'LP_opt', and 'PCL_opt'. 
%
% 
% Primary Contributor: Saeed Azad, PhD

clear; clc; close all;


% Establish the current directory
pdata = EstablishDirectory;

CaseStudies   = {'PCL_opt'};
WaveType      = {'Irregular'};
Regions       = 'EastCoast';
Hydro         = 'MS';
Ctrl          = 'Farm'; 
IP            = 15;
Location_Flag = {'NEAR-SYMMETRICAL'};
Solver        = {'fmincon'};
PowerLimit    = Inf;


pdata.General.n_wec      = 5;

for ii = 1:length(CaseStudies)

    pdata.General.CaseStudy  = CaseStudies{1,ii};
    pdata = CreateDirectories(pdata);
    pdata = ProblemOptions(pdata,Regions,Hydro,Ctrl,IP, Location_Flag{1,1},...
        Solver{1,ii},PowerLimit, WaveType{1,1});

    pdata.Opt.MaxFuncEval    = 1;

    pdata = RunCases(pdata);

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

% Directory path
pdata.General.pathstr  = fileparts(pdata.General.pathstr);

% Go to the current directory
cd(pdata.General.pathstr)

% Add the entire directory to MATLAB's path
addpath(genpath(pdata.General.pathstr))

end