% Run different test cases to ensure they run smoothly with no error
%
%
% Saeed Azad, PhD

clear; clc; close all;

% Establish the current directory
pdata = EstablishDirectory;


CaseStudies   = {'P_opt','C_opt','L_opt','LP_opt','PC_opt','PCL_opt'};
WaveType      = {'Irregular'};
Regions       = 'WestCoast';
Hydro         = 'MS';
Ctrl          = 'Farm'; 
IP            = 15;
Location_Flag = {'NEAR-SYMMETRICAL'};
Solver        = {'surrogateopt','surrogateopt','surrogateopt','surrogateopt','surrogateopt','surrogateopt'};
PowerLimit    = Inf;


pdata.General.n_wec      = 5;

for ii = 1:length(CaseStudies)

    for jj = 1:length(WaveType)

        pdata.General.CaseStudy  = CaseStudies{1,ii};
        pdata = CreateDirectories(pdata);
        pdata = ProblemOptions(pdata,Regions,Hydro,Ctrl,IP, Location_Flag{1,1},...
            Solver{1,ii},PowerLimit, WaveType{1,jj});

        switch upper(pdata.Opt.Solver)
            
            case 'FMINCON'

                pdata.Opt.MaxFuncEval    = 1;

            case 'GA'

                pdata.Opt.PS = 2;
                pdata.Opt.MG = 2; 

            case 'SURROGATEOPT'

                 pdata.Opt.MaxFuncEval = 1;

        end


        pdata = RunCases(pdata);

        fprintf('Case Study %s with wave type %s passed successfully.\n',pdata.General.CaseStudy,WaveType{1,jj})

    end

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