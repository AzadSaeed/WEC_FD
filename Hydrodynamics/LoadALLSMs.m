function pdata = LoadALLSMs(pdata, flag)

% This function loads all the surrogate models to estimate the hydrodynamic
% coefficients. Ten surrogate models were trained using query by committee 
% for each quantity of interest. To enable a robust solution, all these
% ten models are loaded and used in the estimation of the hydrodynamic 
% coefficients. This, however, can be memory-intensive. 
%
%
% Primary Contributor: Saeed Azad, PhD

switch flag

    case 'Non_CMT'

        % Load all 1-WEC models
        ModelSet   = {'A','B','Fe_r', 'Fe_i'};
        Model_dir  = strcat(pdata.General.pathstr, filesep,'SurrogateModels');
        for i = 1:4
            Name = strcat(Model_dir,filesep,'SM_1WEC_',ModelSet{1,i},'_SM_1WEC_lhs_rd_ratio');
            load(Name)
            pdata.SM.Net_1WEC{1,i} = net;
        end

        % Load all 2-WEC models
        ModelSet   = {'A_11','A_12','B_11','B_12','Fe_r_11', 'Fe_i_11'};

        for i = 1:6
            Name = strcat(Model_dir,filesep,'SM_2WEC_',ModelSet{1,i},'_SM_2WEC_Grid_rd_ratio');
            load(Name)
            pdata.SM.Net_2WEC{1,i} = net;
        end

    case 'CMT'
        
        Model_dir  = strcat(pdata.General.pathstr, filesep,'SurrogateModels');
        
        % Load all 1-WEC models
        ModelSet   = {'A','B','Fe_r', 'Fe_i'};

        for i = 1:4
            Name = strcat(Model_dir,filesep,'Models_unscaledInput',filesep,'CMT_1WEC_QBC_',ModelSet{1,i});
            load(Name)
            pdata.SM.Net_1WEC{1,i} = Committee;
        end

        % Load all 2-WEC models
        ModelSet   = {'A_11','A_12','B_11','B_12','Fe_r_11', 'Fe_i_11'};

        for i = 1:6
            Name = strcat(Model_dir,filesep,'Models_unscaledInput',filesep,'CMT_2WEC_QBC_',ModelSet{1,i});
            load(Name)
            pdata.SM.Net_2WEC{1,i} = Committee;
        end

end

end