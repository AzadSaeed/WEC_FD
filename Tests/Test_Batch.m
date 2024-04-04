% Run different test cases to ensure they run smoothly with no error
%
%
% Saeed Azad, PhD

clear; clc; close all;

% Establish the current directory
pdata = EstablishDirectory;

 
Studies       = {'P_opt','C_opt','L_opt','LP_opt','PC_opt','PCL_opt'};
WaveType      = {'Irregular'};
Regions       = {'WestCoast', 'EastCoast'};
Hydro         = {'MS','SM_MBE'};
Ctrl          = {'Farm','Individual'}; 
IP            = {15,10};
Location_Flag = {'NEAR-SYMMETRICAL', 'Far-Symmetrical'}; % For non-layout opt
Solver        = {'GA'};
PowerLimit    = {Inf, 1000};


pdata.General.n_wec      = 5;


for ii = 1:length(Regions)
    for jj = 1:length(Studies)
        for kk = 1:length(Hydro)
            for ll = 1:length(Ctrl)
                for mm = 1:length(IP)
                    for nn = 1:length(Location_Flag)
                        for oo = 1:length(Solver)
                            for pp = 1:length(PowerLimit)
                                for qq = 1:length(WaveType)

                                % Pass data to pdata structure
                                pdata.General.CaseStudy  = Studies{1,jj};

                                % Create directories to store the results
                                pdata = CreateDirectories(pdata);

                                % Precribe all problem data and problem settings
                                pdata = ProblemOptions(pdata,Regions{1,ii},...
                                    Hydro{1,kk},Ctrl{1,ll},IP{1,mm},...
                                    Location_Flag{1,nn}, Solver{1,oo},... ...
                                    PowerLimit{1,pp},WaveType{1,qq});

                                switch pdata.Opt.Solver

                                    case 'fmincon'

                                        pdata.Opt.MaxFuncEval    = 1;

                                    case 'GA'

                                        pdata.Opt.PS = 2;
                                        pdata.Opt.MG = 1;

                                    case 'SurrogateOpt'

                                        pdata.Opt.MaxFuncEval = 3;

                                end


                                pdata = RunCases(pdata);

                                close all;

                                % To avoid having multiple interactive
                                % sessions
                                % delete(gcp('nocreate'))

                                if strcmpi(Studies{1,jj},'P_opt') || strcmpi(Studies{1,jj},'C_opt') || strcmpi(Studies{1,jj},'PC_opt')
                                    fprintf(['Batch Case Study %s with wave type %s in region %s, \n'...
                                        'using %s for hydrodynamics, with %s control, %d starting point,\n'...
                                        'solver %s, power limit %d, at layout %s passed successfully.\n'], ...
                                        Studies{1,jj},WaveType{1,qq},Regions{1,ii},Hydro{1,kk},Ctrl{1,ll},...
                                        IP{1,mm},Solver{1,oo}, PowerLimit{1,pp},Location_Flag{1,nn})
                                elseif strcmpi(Studies{1,jj},'L_opt') || strcmpi(Studies{1,jj},'LP_opt') || strcmpi(Studies{1,jj},'PCL_opt')
                                    fprintf(['Batch Case Study %s with wave type %s in region %s, \n'...
                                        'using %s for hydrodynamics, with %s control, %d starting point, \n'...
                                        'and solver %s, power limit %d passed successfully.\n'],...
                                        Studies{1,jj},WaveType{1,qq},Regions{1,ii},Hydro{1,kk},...
                                        Ctrl{1,ll},IP{1,mm},Solver{1,oo}, PowerLimit{1,pp})
                                end

                                end
                            end
                        end
                    end
                end
            end
        end
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