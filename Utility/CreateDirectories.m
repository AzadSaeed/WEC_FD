function pdata = CreateDirectories(pdata)
% This function creates directories for various case studies. To avoid
% overwriting the solutions, the directories are created based on the
% following criteria in order:
%
% ---- Solution
% ------------- Case Study
% ------------------------ Number of WECs
% --------------------------------------- Wave Type
% -------------------------------------------------- results or Plots
%
%
% Contributor: Saeec Azad, PhD


%-------------------------------------------------------------------------%
% Working Directory
pdata.General.working_dir = pwd;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Storage Directory
pdata.General.Storage_dir = strcat('Solution',filesep,pdata.General.CaseStudy);
if ~isfolder(pdata.General.Storage_dir)
    mkdir(pdata.General.Storage_dir)
end
%-------------------------------------------------------------------------%


pdata = CreateDir(pdata);



end



function pdata = CreateDir(pdata)

n_wec = pdata.General.n_wec;

if ~strcmpi(pdata.General.CaseStudy,'SIM')

    %-------------------------------------------------------------------------%
    % Directory based on the number of WECs
    pdata.General.nwec_dir = strcat(pdata.General.Storage_dir, filesep, string(n_wec),'_WEC');
    if ~isfolder(pdata.General.nwec_dir)
        mkdir(pdata.General.nwec_dir)
    end
    %-------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % Results Directory for simulation with FD and regular waves
    pdata.General.paths.Sol_R = strcat(pdata.General.nwec_dir, filesep, 'Regular');
    if ~isfolder(pdata.General.paths.Sol_R)
        mkdir(pdata.General.paths.Sol_R)
    end
    %-------------------------------------------------------------------------%

    % Create separate plot directories for better organization
    %-------------------------------------------------------------------------%
    % Plot Directory
    pdata.General.plot_dir_R = strcat(pdata.General.paths.Sol_R,filesep,'Plots');
    if ~isfolder(pdata.General.plot_dir_R)
        mkdir(pdata.General.plot_dir_R)
    end
    %-------------------------------------------------------------------------%

    % Create separate plot directories for better organization
    %-------------------------------------------------------------------------%
    % Plot Directory
    pdata.General.sol_dir_R = strcat(pdata.General.paths.Sol_R,filesep,'Results');
    if ~isfolder(pdata.General.sol_dir_R)
        mkdir(pdata.General.sol_dir_R)
    end
    %-------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % Results Directory for simulation with FD and irregular waves
    pdata.General.paths.Sol_IR = strcat(pdata.General.nwec_dir, filesep, 'Irregular');
    if ~isfolder(pdata.General.paths.Sol_IR)
        mkdir(pdata.General.paths.Sol_IR)
    end
    %-------------------------------------------------------------------------%


    %-------------------------------------------------------------------------%
    % Plot Directory
    pdata.General.plot_dir_IR  = strcat(pdata.General.paths.Sol_IR,filesep,'Plots');
    if ~isfolder(pdata.General.plot_dir_IR)
        mkdir(pdata.General.plot_dir_IR)
    end
    %-------------------------------------------------------------------------%

    %-------------------------------------------------------------------------%
    % Plot Directory
    pdata.General.sol_dir_IR = strcat(pdata.General.paths.Sol_IR,filesep,'Results');
    if ~isfolder(pdata.General.sol_dir_IR)
        mkdir(pdata.General.sol_dir_IR)
    end
    %-------------------------------------------------------------------------%

    % if strcmpi(pdata.General.CaseStudy,'PCL_opt')
    % 
    %     % Nested
    %     pdata.General.sol_dir_IR_OptArch_N = strcat(pdata.General.sol_dir_IR,filesep,'Nested');
    %     if ~isfolder(pdata.General.sol_dir_IR_OptArch_N)
    %         mkdir(pdata.General.sol_dir_IR_OptArch_N)
    %     end
    % 
    %     % Simultaneous
    %     pdata.General.sol_dir_IR_OptArch_S = strcat(pdata.General.sol_dir_IR,filesep,'Simultaneous');
    %     if ~isfolder(pdata.General.sol_dir_IR_OptArch_S)
    %         mkdir(pdata.General.sol_dir_IR_OptArch_S)
    %     end
    % 
    % 
    % end

elseif strcmpi(pdata.General.CaseStudy,'SIM')

    pdata.General.Solution_dir = strcat(pdata.General.Storage_dir,filesep,pdata.SimCaseStudy);
    if ~isfolder(pdata.General.Solution_dir)
        mkdir(pdata.General.Solution_dir)
    end
    pdata.General.plot_dir = pdata.General.Solution_dir;


end
end
