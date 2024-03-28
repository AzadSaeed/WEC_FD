function ACC_Runs(pdata)
% This function runs the cases investigated in "Site-dependent Solutions of
% Wave Energy Converter Farms with Surrogate Models, Control Co-design, and
% Layout Optimization" accepted in American Control Conference 2024 by
% Saeed Azad, Daniel R. Herber, Suraj Khanal, Gaofeng Jia
%
% To run, simply put this function under the BatchStudies case in the
% Main_WEC_FD file.
% 
%
% Contributor: Saeed Azad, PhD


ACC_CaseStudy = 'CaseStudyIII';


switch upper(ACC_CaseStudy)

    case 'CASESTUDYI'

        %-----------------------------------------------------------------%
        %------------------------- Case Study I --------------------------%
        %-----------------------------------------------------------------%

        Studies       = {'LP_opt'};
        Regions       = {'AlaskaCoasts','EastCoast','PacificIslands','WestCoast'};
        Hydro         = {'SM_MBE'};
        IP            = {15};
        Location_Flag = {'NEAR-SYMMETRICAL'};
        Solver        = {'GA'};
        PowerLimit    = {Inf};
        Ctrl          = {'Farm'};
        WaveType      = {'Irregular'};

        % Run Batch Studies
        pdata = RunBatchStudies(Regions,Studies,Hydro,IP,Ctrl,Location_Flag,...
            Solver,PowerLimit, WaveType, pdata);


    case 'CASESTUDYII'

        %-----------------------------------------------------------------%
        %----------------------- Case Study II ---------------------------%
        %-----------------------------------------------------------------%

        Regions       = {'AlaskaCoasts','EastCoast','PacificIslands','WestCoast'};
        Studies       = {'PCL_opt'};
        Hydro         = {'SM_MBE'};
        IP            = {15};
        Location_Flag = {'NEAR-SYMMETRICAL'};
        Solver        = {'GA'};
        PowerLimit    = {Inf};
        Ctrl          = {'Farm'};
        WaveType      = {'Irregular'};

        % Run batch studies
        pdata = RunBatchStudies(Regions,Studies,Hydro,IP,Ctrl,Location_Flag,...
            Solver, PowerLimit,WaveType, pdata);


    case 'CASESTUDYIII'
        
        %-----------------------------------------------------------------%
        %----------------------- Case Study III --------------------------%
        %-----------------------------------------------------------------%

        Regions       = {'AlaskaCoasts','EastCoast','PacificIslands','WestCoast'};
        Studies       = {'PCL_opt'};
        Hydro         = {'SM_MBE'};
        IP            = {15};
        Location_Flag = {'NEAR-SYMMETRICAL'};
        Solver        = {'GA'};
        PowerLimit    = {Inf};
        Ctrl          = {'Individual'};
        WaveType      = {'Irregular'};

        % Run batch studies
        pdata = RunBatchStudies(Regions,Studies,Hydro,IP,Ctrl, Location_Flag,...
            Solver, PowerLimit, WaveType, pdata);

end
end