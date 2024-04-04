% This script:
% 
% 1. loads the solutions from Case II and Case III associated with 
% ACC paper ('Site-dependent Solutions of Wave Energy Converter Farms with
% Surrogate Models, Control Co-design, and Layout Optimization') 
%
% 2. Creates a large number of randomly selected layouts and runs the
% simulations with fixed optimized plant and control 
%
%
% 3. Analyzes how the performance of the farm compares with randomly 
% selected samples 
%
% 
% Primary Contributor: Saeed Azad. PhD


clear; clc; close all;


% Establish the current directory
pdata = EstablishDirectory;

% Case Study to run/load
% CaseStudies = {'CaseStudyI','CaseStudyII','CaseStudyIII'};
CaseStudies = {'CaseStudyII'};
Locations = {'AlaskaCoasts','EastCoast','PacificIslands','WestCoast'};
Method = 'Simulations';


for ii = 1:length(CaseStudies)
    for jj = 1:length(Locations)


        CaseStudy = CaseStudies{1,ii};
        Loc = Locations{1,jj};
        
        % Load ACC Solution
        sol = LoadACCSolutions(CaseStudy,Loc, pdata);


        % Prescribe the required specifications
        pdata.General.n_wec     = 5;
        Regions                 = Locations{1,jj};
        pdata.HydroFlag         = 'MS';
        PowerLimit              = Inf;
        pdata.General.CaseStudy = 'Sim';


        % Create directories to store the results
        % pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,Regions,pdata.HydroFlag,...
            sol.pdata.Opt.CtrlFlag,15,[],sol.pdata.Opt.Algorithm,PowerLimit,...
            'Irregular');

        % Load all surrogate models once-if needed
        if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
            % 'CMT' or 'Non-CMT'
            pdata = LoadALLSMs(pdata,'CMT');
        end
       

        % Update directories 
        pdata.General.Solution_dir = strcat(pdata.General.pathstr,filesep,'Solution',filesep, 'ACC', filesep, CaseStudy, filesep, 'Analysis');
        pdata.General.plot_dir = strcat(pdata.General.pathstr,filesep,'Solution',filesep, 'ACC', filesep, CaseStudy, filesep, 'Analysis');
        if ~isfolder(pdata.General.Solution_dir)
            mkdir(pdata.General.Solution_dir)
        end
        if ~isfolder(pdata.General.plot_dir)
            mkdir(pdata.General.plot_dir)
        end


        switch upper(Method)

            case 'SIMULATIONS'

                % Farm area to investigate 
                farm_limits_x = [sol.pdata.Opt.lbl(1,1), sol.pdata.Opt.ubl(1,1)];
                farm_limits_y = [sol.pdata.Opt.lbl(2,1), sol.pdata.Opt.ubl(2,1)];

                % Number of random layouts to generate
                N_l = 2500;

                % Prescribe Radius
                pdata.WEC.Radius  = sol.pdata.Results.Radius;
                pdata.WEC.RD      = sol.pdata.Results.RD;
                pdata.WEC.Draft   = pdata.WEC.Radius./pdata.WEC.RD;

                % Prescribe PTO parameters
                pdata.WEC.PTO_damping    = sol.pdata.Results.PTO_damping;
                pdata.WEC.PTO_stiffness  = sol.pdata.Results.PTO_stiffness;

                % Randomly create layouts
                s0Lx = 2651; s0Ly = 12584;
                x_ = NaN(pdata.General.n_wec-1,N_l);
                y_ = NaN(pdata.General.n_wec-1,N_l);
                kk = 0;

                for i = 1:N_l

                    % Make it reproducible
                    rng(s0Lx+kk)
                    
                    x_(:,i) = farm_limits_x(1,1) + ...
                        (farm_limits_x(1,2)-farm_limits_x(1,1))*rand([1,pdata.General.n_wec-1]);
                    
                    % Make it reproducible
                    rng(s0Ly+kk)

                    y_(:,i) = farm_limits_y(1,1) + ...
                        (farm_limits_y(1,2)-farm_limits_y(1,1) )*rand([1,pdata.General.n_wec-1]);
                    kk = kk+1;
                end

                % Assign WEC locations with first WEC at the center
                pdata.WEC.x = [zeros(1,N_l);x_];
                pdata.WEC.y = [zeros(1,N_l);y_];

                Dist(1,1)   = 0;
                Angle(1,1)  = 0;

                % Define the reference WEC
                ref     = [0,0];

                for j = 1:N_l
                    xl = pdata.WEC.x(:,j);
                    yl = pdata.WEC.y(:,j);
                    for i=2:pdata.General.n_wec
                        delta_x         = xl(i,1) - ref(1,1);
                        delta_y         = yl(i,1) - ref(1,2);
                        Dist(i,j)       = sqrt(delta_x.^2 + delta_y.^2);
                        Angle(i,j)      = atan2(delta_y,delta_x);
                    end
                end

                pdata.WEC.Distance      = Dist;
                pdata.WEC.Angle         = Angle;

                % Check to ensure non-overlapping layout
                pdata.N_l = N_l;
                pdata = CheckSafetyDis(pdata);


                % Update number of simulations
                pdata.N_l = length(pdata.WEC.x(1,:));

                % Run simulations

                X   = pdata.WEC.x;
                Y   = pdata.WEC.y;
                Dis = pdata.WEC.Distance;
                Ang = pdata.WEC.Angle;

                % Dummy variable - needed for parallel implementation
                P_DATA = pdata;

                % number of simulations
                N_l  = pdata.N_l;


                % Pre-assign
                Power = NaN(1,N_l);
                Q_fac  = cell(1,N_l);
                q_fac  = cell(1,N_l);

                t1 = tic;

                % Perform simulation
                nw = 3;
                parfor (i = 1:N_l,nw)

                    % Broadcast variable
                    P_DATA               = pdata;

                    % Layout variables
                    P_DATA.WEC.x         = X(:,i);
                    P_DATA.WEC.y         = Y(:,i);
                    P_DATA.WEC.Distance  = Dis(:,i);
                    P_DATA.WEC.Angle     = Ang(:,i);


                    % run simulation
                    Sim_sol         = run_simulation(P_DATA);

                    % Power [W]
                    Power(1,i)     = sum(sum(Sim_sol.Results.Power.P_total_vec));


                    % Calculate q-factor
                    [Q_factor, q_factor]     = calculate_q_factor(Sim_sol);

                    Q_fac{1,i} = Q_factor;
                    q_fac{:,i} = q_factor;

                end

                Q.Q_fac = Q_fac;
                Q.q_fac = q_fac;

                Volume = pi*(pdata.WEC.Radius.^2)*2.*pdata.WEC.Draft;
                Obj    = Power/Volume; % W/m^3

                t2 = toc(t1);

                pdata.Results.Q = Q;
                pdata.Results.Power = Power;
                pdata.Results.Obj = Obj;

                savename = strcat(pdata.General.Solution_dir,filesep,'Analysis_',Loc);
                save(savename,'pdata')


            case 'PERTURB'
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
pdata.General.pathstr  = fileparts(pdata.General.pathstr);

% Go to the current directory
cd(pdata.General.pathstr)

% Add the entire directory to MATLAB's path
addpath(genpath(pdata.General.pathstr))

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function sol = LoadACCSolutions(CaseStudy,Loc, pdata)

% ACC Solution directory
switch upper(CaseStudy)

    case 'CASESTUDYI'


        sol_dir = strcat(pdata.General.pathstr,filesep, 'Solution', filesep,...
            'ACC', filesep, CaseStudy,filesep, 'SM_MBE', filesep,Loc, ...
            filesep,'Solution_GA_LP_opt_500_SM_MBE_rng_15_Plim_None');

    case 'CASESTUDYII'

        sol_dir = strcat(pdata.General.pathstr,filesep, 'Solution', filesep,...
            'ACC', filesep, CaseStudy,filesep, Loc, ...
            filesep,'Solution_GA_PCL_opt_50_SM_MBE_Farm_rng_15_Plim_None');

    case 'CASESTUDYIII'

        sol_dir = strcat(pdata.General.pathstr,filesep, 'Solution', filesep,...
            'ACC', filesep, CaseStudy,filesep, Loc, ...
            filesep,'Solution_GA_PCL_opt_50_SM_MBE_Individual_rng_15_Plim_None');
end


% Load soluiton
sol = load(sol_dir);

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = CheckSafetyDis(pdata)

idx_ = 1:pdata.General.n_wec;
C    = nchoosek(idx_,uint16(2));
idr_vec = [];
id_c = 1;

for i = 1:pdata.N_l
    
    xx  = pdata.WEC.x(:,i);
    yy  = pdata.WEC.y(:,i);
    xl  = [xx,yy];
    rr  = pdata.WEC.Radius;


    for ii = 1:length(C(:,1))

        p_body  = C(ii,1);
        q_body  = C(ii,2);
        x_p     = xl(p_body,1);
        y_p     = xl(p_body,2);
        x_q     = xl(q_body,1);
        y_q     = xl(q_body,2);
        l_pq    = sqrt((x_p - x_q)^2 + (y_p - y_q)^2);
        c       = 2*rr + pdata.General.SafeDist - l_pq;

        if c > 0
           idr_vec(id_c) = i;
           id_c = id_c + 1;
        end

    end

end


if ~isempty(idr_vec)

    pdata.WEC.x(:,idr_vec) = [];
    pdata.WEC.y(:,idr_vec) = [];

    pdata.WEC.Distance(:,idr_vec)  = [];
    pdata.WEC.Angle(:,idr_vec)     = [];

end


end
