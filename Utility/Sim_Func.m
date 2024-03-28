function pdata = Sim_Func(pdata)
% This function defines the required conditions for various simulations.
% Some pre-defined simulations are described below.
%
%
% Case 1_Regular: 1 WEC with fixed radius and draft, at the center of the
%           coordinate system with prescribed values for PTO parameters
%
%
% Case_2_Regular: 2 WECs with fixed radius and draft, far from eachother (distance 
%           900~m and angle pi/6), with prescribed PTO parameters
%
%
% Case_1_Irregular: 2 WECs with fixed radius and draft, clsoe to eachother (distance 
%           20~m and angle pi/6), with prescribed PTO parameters
%
%
% Case_2_Irregular: 2 WECs with fixed radius and draft, arbitrary array, with 
%           prescribed PTO parameters
%
%
% Case_Farm_landscape 2 WECs at a different grid locations with prescribed PTO
%           parameters to create a surf plot of power and q-factor
%
%
% Case_layout_sample: case 'batch_simulations': 5 WECS located at random with 
%           randomized PTO and plant for each case. The goal is to understand
%           the error in terms of power from the surrogate models
%           case 'run_optimized':Run simulation from optimized results 
%
%
% Case_Radius_v_kPTO: runs a series of simulations with a fixed slenderness
%           ratio, and various, PTO parameters and radii
%
%
%
% Case_Load_Optim: Run simulation from optimized results for various cases
%
%
% Primary Contributer: Saeed Azad, PhD




pdata.SimCaseStudy = 'Case_Load_Optim';

switch pdata.SimCaseStudy

    case 'Case_1_Regular'

        pdata.General.n_wec     = 1;
        pdata.HydroFlag         = 'MS';
        PowerLimit              = Inf;
        pdata.Opt.ObjForm       = 'Power';

        % Create directories to store the results
         pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,[],pdata.HydroFlag,[],[],[],[],PowerLimit, 'Regular');

        pdata.WEC.Radius        = 4.24;
        pdata.WEC.RD            = 8.48;
        pdata.WEC.Draft         = pdata.WEC.Radius /pdata.WEC.RD;
        pdata.WEC.Distance      = 0;
        pdata.WEC.Angle         = 0;
        pdata.WEC.x             = (pdata.WEC.Distance.*cos(pdata.WEC.Angle));
        pdata.WEC.y             = (pdata.WEC.Distance.*sin(pdata.WEC.Angle));
        pdata.WEC.PTO_damping   = 5*10^5*ones(1,pdata.General.n_wec);
        pdata.WEC.PTO_stiffness = -5*10^3*ones(1,pdata.General.n_wec);
        pdata.General.FindOptimalPTO = 0;
        pdata.Opt.SatFlag       = 0;

    case 'Case_2_Regular'


        pdata.General.n_wec = 1; 

        switch pdata.General.n_wec

            case 1

                % Assign WEC locations with first WEC at the center
                pdata.WEC.x = 0;
                pdata.WEC.y = 0;

            case 2

                % Assign WEC locations with first WEC at the center
                pdata.WEC.x = [0;0];
                pdata.WEC.y = [0;8];

            case 3

                % Assign WEC locations with first WEC at the center
                pdata.WEC.x = [0;3;1.5];
                pdata.WEC.y = [0;0;7.5];


        end

        pdata.HydroFlag         = 'MS';
        PowerLimit              = Inf;
        pdata.Opt.ObjForm       = 'Power';

        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,[],pdata.HydroFlag,[],[],[],[],PowerLimit, 'Regular');

        pdata.WEC.Radius          = 1;
        pdata.WEC.RD              = 1;
        pdata.WEC.Draft           = pdata.WEC.Radius /pdata.WEC.RD;
        pdata.General.depth       = 8;
        k                         = 0.4/pdata.WEC.Radius;
        omega                     = sqrt(pdata.General.g*k*tanh(k*pdata.General.depth));
        Tp                        = (2*pi)/omega;
        pdata.WaveData.Te_regular = Tp;

        A = 1;
        H = 2*A;
        Hs = sqrt(2)*H;
        pdata.WaveData.Hs_regular = Hs;

        pdata  = CallWaveSpectrum(pdata);

        pdata.WEC.PTO_stiffness = zeros(1,pdata.General.n_wec);
        pdata.WEC.PTO_damping   = 5*10^5*ones(1,pdata.General.n_wec);

        % pdata.Results.PTO_stiffness = pdata.WEC.PTO_stiffness;
        % pdata.Results.Radius        = pdata.WEC.Radius;
        % pdata.Results.Draft         = pdata.WEC.Draft;
        % pdata.Results.x             = pdata.WEC.x;
        % pdata.Results.y             = pdata.WEC.y;
        pdata.General.FindOptimalPTO = 0;
        pdata.Opt.SatFlag            = 0;

    case 'Case_1_Irregular'

        pdata.General.n_wec     = 1;
        Regions                 = 'WestCoast';
        pdata.HydroFlag         = 'SM_MBE';
        PowerLimit              = Inf;

        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,Regions,pdata.HydroFlag,[],[],[],[],PowerLimit, 'Irregular');

        % Load all surrogate models once-if needed
        if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
            % 'CMT' or 'Non-CMT'
            pdata = LoadALLSMs(pdata,'CMT');
        end


        pdata.WEC.Radius        = 10;
        pdata.WEC.RD            = 2;
        pdata.WEC.Draft         = pdata.WEC.Radius /pdata.WEC.RD;
        pdata.WEC.Distance      = 0;
        pdata.WEC.Angle         = 0;
        pdata.WEC.x             = (pdata.WEC.Distance.*cos(pdata.WEC.Angle));
        pdata.WEC.y             = (pdata.WEC.Distance.*sin(pdata.WEC.Angle));
        pdata.WEC.PTO_damping   = 8461.53189662730;
        pdata.WEC.PTO_stiffness = -366990.727560625;

    case 'Case_2_Irregular'

        pdata.General.n_wec     = 2;
        Regions                 = 'WestCoast';
        pdata.HydroFlag         = 'SM_MBE';
        PowerLimit              = Inf;

        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,Regions,pdata.HydroFlag,[],[],[],[],PowerLimit, 'Irregular');

        % Load all surrogate models once-if needed
        if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
            % 'CMT' or 'Non-CMT'
            pdata = LoadALLSMs(pdata,'CMT');
        end

        pdata.WEC.Radius        = 5;
        pdata.WEC.RD            = 10;
        pdata.WEC.Draft         = pdata.WEC.Radius /pdata.WEC.RD;
        % pdata.WEC.Distance      = [0;2*pdata.WEC.Radius+pdata.General.SafeDist];
        pdata.WEC.Distance      = [0;900];
        pdata.WEC.Angle         = [0;pi/6];
        pdata.WEC.x             = (pdata.WEC.Distance.*cos(pdata.WEC.Angle));
        pdata.WEC.y             = (pdata.WEC.Distance.*sin(pdata.WEC.Angle));
        pdata.WEC.PTO_damping   = [1, 1];
        pdata.WEC.PTO_stiffness = [-10^6, -10^6];

    case 'Case_3_Irregular'

        pdata.General.n_wec     = 5;
        Regions                 = 'WestCoast';
        pdata.HydroFlag         = 'SM_MBE';
        PowerLimit              = Inf;

        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,Regions,pdata.HydroFlag,[],[],[],[],PowerLimit, 'Irregular');

        % Load all surrogate models once-if needed
        if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
            % 'CMT' or 'Non-CMT'
            pdata = LoadALLSMs(pdata,'CMT');
        end

        pdata.WEC.Radius        = 3.59947702216808;
        pdata.WEC.RD            = 6.95077708704819;
        pdata.WEC.Draft         = pdata.WEC.Radius /pdata.WEC.RD;

        % WEC location
        pdata.WEC.x = [0; 143.018090425836; 102.900204606563; 63.3050138064852; 62.8852056020629];
        pdata.WEC.y = [0; 21.2706134758765; -42.3890646968897; -117.607185266508; 103.417784001519];

        pdata.WEC.PTO_damping   = 8461.53189662730*ones(1,pdata.General.n_wec);
        pdata.WEC.PTO_stiffness = -366990.727560625*ones(1,pdata.General.n_wec);
        pdata.Opt.SatFlag       = 0;
        pdata.General.FindOptimalPTO = 0;


    case 'Case_Farm_landscape'

        pdata.General.n_wec     = 2;
        Regions                 = 'WestCoast';
        pdata.HydroFlag         = 'SM_MBE';
        PowerLimit              = Inf;

        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,Regions,pdata.HydroFlag,[],[],[],[],PowerLimit, 'Irregular');

        % Load all surrogate models once-if needed
        if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
            % 'CMT' or 'Non-CMT'
            pdata = LoadALLSMs(pdata,'CMT');
        end

        n1 = 2;
        n2 = 3;
        pdata.WEC.Radius        = 5;
        pdata.WEC.RD            = 10;
        pdata.WEC.Draft         = pdata.WEC.Radius /pdata.WEC.RD;
        pdata.WEC.x_vec_p       = logspace(log10(6),log10(300),n1);
        pdata.WEC.y_vec_p       = logspace(log10(6),log10(300),n2);
        pdata.WEC.x_vec_n       = -pdata.WEC.x_vec_p;
        pdata.WEC.y_vec_n       = -pdata.WEC.y_vec_p;
        x_vec                   = [fliplr(pdata.WEC.x_vec_n),pdata.WEC.x_vec_p];
        y_vec                   = [flip(pdata.WEC.y_vec_n),pdata.WEC.y_vec_p];
        [x_mat,y_mat]           = meshgrid(x_vec,y_vec);
        X_vec                   = x_mat(:);
        Y_vec                   = y_mat(:);
        pdata.WEC.x             = [zeros(length(X_vec),1),X_vec]';
        pdata.WEC.y             = [zeros(length(X_vec),1),Y_vec]';
        pdata.WEC.Distance      = sqrt(pdata.WEC.x.^2 + pdata.WEC.y.^2);
        pdata.WEC.Angle         = atan2(pdata.WEC.y,pdata.WEC.x);
        pdata.WEC.PTO_damping   = [1000, 1000];
        pdata.WEC.PTO_stiffness = [-10^5, -10^5];
        pdata.n_sim             = length(x_vec)*length(y_vec); % number of simulations
        pdata.General.plotflag = 0;

    case 'Case_layout_sample'

        % Cases to run: 'batch_simulations' or 'run_optimized'
        pdata.General.run_case = 'run_optimized';
        pdata.General.n_wec     = 5;

        Regions                 = 'WestCoast';
        pdata.HydroFlag         = 'SM_MBE';
        PowerLimit              = Inf;

        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,Regions,pdata.HydroFlag,[],[],[],[],PowerLimit, 'Irregular');
        pdata.Opt.CtrlFlag  = 'Farm';

        % Load all surrogate models once-if needed
        if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
            % 'CMT' or 'Non-CMT'
            pdata = LoadALLSMs(pdata,'CMT');
        end

         pdata.General.plotflag = 0;


        switch pdata.General.run_case

            case 'batch_simulations'

                % Number of layouts to investigate between [-150 150] in x & y
                farm_limits = [-150, 150];
                N_l = 5;

                pdata.Opt.R_min          = 0.5;
                pdata.Opt.R_max          = 10;

                % Radius to draft ratio  [m]
                pdata.Opt.RD_min         = 0.2;
                pdata.Opt.RD_max         = 10;

                % Draft [m] -- defined to be half of the overall height
                pdata.Opt.D_min          = 0.5;
                pdata.Opt.D_max          = 20;

                % PTO [Damping (N.s/m), stiffness (N/m)]
                pdata.Opt.k_PTO_min      = -5*10^5;
                pdata.Opt.k_PTO_max      = 5*10^5;
                pdata.Opt.B_PTO_min      = 0;
                pdata.Opt.B_PTO_max      = 5*10^5;

                % Plant design - selected randomly
                rng(2154);
                r_ = pdata.Opt.R_min + (pdata.Opt.R_max-pdata.Opt.R_min )*rand([1,N_l]);

                rng(48511);
                ar_ = pdata.Opt.RD_min + (pdata.Opt.RD_max-pdata.Opt.RD_min)*rand([1,N_l]);


                pdata.WEC.Radius        = r_;
                pdata.WEC.RD            = ar_;
                pdata.WEC.Draft         = pdata.WEC.Radius./pdata.WEC.RD;

                % Layout design - selected randomly
                s0Lx = 2651; s0Ly = 12584;
                x_ = NaN(pdata.General.n_wec-1,N_l);
                y_ = NaN(pdata.General.n_wec-1,N_l);
                kk = 0;

                for i = 1:N_l
                    rng(s0Lx+kk)
                    x_(:,i) = farm_limits(1,1) + (farm_limits(1,2)-farm_limits(1,1) )*rand([1,pdata.General.n_wec-1]);
                    rng(s0Ly+kk)
                    y_(:,i) = farm_limits(1,1) + (farm_limits(1,2)-farm_limits(1,1) )*rand([1,pdata.General.n_wec-1]);
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


                % control design - selected randomly
                if strcmpi(pdata.Opt.CtrlFlag,'Individual')
                    s0Lb = 875215; s0Lk = 95812;
                    B_PTO = NaN(pdata.General.n_wec,N_l);
                    K_PTO = NaN(pdata.General.n_wec,N_l);
                    kk = 0;
                    for i = 1:N_l
                        rng(s0Lb+kk)
                        B_PTO(:,i) = (pdata.Opt.B_PTO_min + (pdata.Opt.B_PTO_max-pdata.Opt.B_PTO_min)*rand([1,pdata.General.n_wec]));
                        rng(s0Lk+kk)
                        K_PTO(:,i) = pdata.Opt.k_PTO_min + (pdata.Opt.k_PTO_max -pdata.Opt.k_PTO_min)*rand([1,pdata.General.n_wec]);
                        kk = kk+1;
                    end
                elseif strcmpi(pdata.Opt.CtrlFlag,'FARM')
                    s0Lb = 875215; s0Lk = 95812;
                    B_PTO = NaN(1,N_l);
                    K_PTO = NaN(1,N_l);
                    kk = 0;
                    for i = 1:N_l
                        rng(s0Lb+kk)
                        B_PTO(:,i) = (pdata.Opt.B_PTO_min + (pdata.Opt.B_PTO_max-pdata.Opt.B_PTO_min)*rand([1,1]));
                        rng(s0Lk+kk)
                        K_PTO(:,i) = pdata.Opt.k_PTO_min + (pdata.Opt.k_PTO_max -pdata.Opt.k_PTO_min)*rand([1,1]);
                        kk = kk+1;
                    end
                    B_PTO = repmat(B_PTO,[pdata.General.n_wec,1]);
                    K_PTO = repmat(K_PTO,[pdata.General.n_wec,1]);
                end


                pdata.WEC.PTO_damping   = B_PTO';
                pdata.WEC.PTO_stiffness = K_PTO';
                pdata.N_l = N_l;


                % Check to ensure non-overlapping layout
                pdata = CheckSafetyDis(pdata);

                % Remove points where Draft is not within limits
                pdata = CheckDraft(pdata);

                % Update number of simulations
                pdata.N_l = length(pdata.WEC.Radius);

                % savename = strcat(pdata.General.Storage_dir, filesep, 'Case_', pdata.SimCaseStudy,filesep, 'case_6_input_data');
                % In = pdata.WEC;
                % save(savename,'In');

            case 'run_optimized'

                try

                    sol = load('xxx');
                    pdata.WEC.Radius        = sol.pdata.Results.Radius;
                    pdata.WEC.RD            = sol.pdata.Results.RD;
                    pdata.WEC.Draft         = sol.pdata.Results.Draft;
                    pdata.WEC.PTO_damping   = sol.pdata.Results.PTO_damping;
                    pdata.WEC.PTO_stiffness = sol.pdata.Results.PTO_stiffness;
                    pdata.WEC.x             = sol.pdata.Results.x;
                    pdata.WEC.y             = sol.pdata.Results.y;
                    pdata.WEC.Distance      = sol.pdata.Results.Distance;
                    pdata.WEC.Angle         = sol.pdata.Results.Angle;
                    pdata.N_l               = 1; % number of simulations

                catch

                    error('Could not load solution for the simulation.')

                end
        end

    case 'Case_Radius_v_kPTO'


        pdata.General.n_wec = 5;
        PowerLimit              = Inf;

        % Location flag and hydrodynamics
        pdata.General.Location_Flag = 'NEAR-SYMMETRICAL';
        pdata.HydroFlag             = 'MS';


        % Create directories to store the results
        pdata = CreateDirectories(pdata);

        % Precribe all problem data and problem settings
        pdata = ProblemOptions(pdata,'WestCoast',pdata.HydroFlag,[],[],[],[],PowerLimit, 'Irregular');
        pdata.Opt.CtrlFlag  = 'Farm';

        % Fixed aspect ratio
        RD      = 1;

        % Number of points in radius dimension
        n1      = 3;
        pdata.General.n1 = n1;

        % Number of points in k_pto dimension
        n2      = 3;
        pdata.General.n2 = n2;

        % Number of points in B_pto dimension
        n3      = 4;
        pdata.General.n3 = n3;

        % Radius bounds
        pdata.Opt.R_min          = 0.5;
        pdata.Opt.R_max          = 10;

        % Radius to draft ratio  [m]
        pdata.Opt.RD_min         = 0.2;
        pdata.Opt.RD_max         = 10;

        % Draft [m] -- defined to be half of the overall height
        pdata.Opt.D_min          = 0.5;
        pdata.Opt.D_max          = 20;

        % PTO [Damping (N.s/m), stiffness (N/m)]
        pdata.Opt.k_PTO_min      = -5*10^5;
        pdata.Opt.k_PTO_max      = 5*10^5;
        pdata.Opt.B_PTO_min      = 0;
        pdata.Opt.B_PTO_max      = 5*10^5;

        % Radius vector
        r_vec = logspace(log10(pdata.Opt.R_min),log10(pdata.Opt.R_max),n1);

        % B_pto vector
        B_vec = linspace(pdata.Opt.B_PTO_min+0.1,pdata.Opt.B_PTO_max,n2);

        % K_pto vector
        K_vec = linspace(pdata.Opt.k_PTO_min,pdata.Opt.k_PTO_max,n3);


        % Combine the two vectors
        rk    = {r_vec, B_vec, K_vec};
        C_1   = cell(size(rk));
        [C_1{:}] = ndgrid(rk{:});
        C_1      = cellfun(@(x) x(:),C_1,'un',0);
        RK       = [C_1{1,1:3}];

        % Specify the locatioo
        [x,y] = Locations(pdata.General.Location_Flag , pdata.General.n_wec);
        pdata.WEC.x = x;
        pdata.WEC.y = y;

        pdata.WEC.Radius        = RK(:,1);
        pdata.WEC.RD            = RD;
        pdata.WEC.Draft         = pdata.WEC.Radius /pdata.WEC.RD;

        pdata.WEC.PTO_damping     = RK(:,2)*ones(1,pdata.General.n_wec);
        pdata.WEC.PTO_stiffness   = RK(:,3)*ones(1,pdata.General.n_wec);

        % number of simulations
        pdata.n_sim             = n1*n2*n3;



    case 'Case_Load_Optim'

        try

            pdata.General.n_wec     = 5;
            Regions                 = 'EastCoast';
            pdata.HydroFlag         = 'MS';
            PowerLimit              = Inf;

            % Create directories to store the results
            pdata = CreateDirectories(pdata);

            % Precribe all problem data and problem settings
            pdata = ProblemOptions(pdata,Regions,pdata.HydroFlag,[],[],[],[],PowerLimit, 'Irregular');
            pdata.Opt.CtrlFlag  = 'Farm';


            % Load all surrogate models once-if needed
            if strcmpi(pdata.HydroFlag,'SM_MBE') || strcmpi(pdata.HydroFlag,'MS_MBE')
                % 'CMT' or 'Non-CMT'
                pdata = LoadALLSMs(pdata,'CMT');
            end

            
            Location = Regions;
            Alg      = 'GA';
            study    = 'PCL_opt';
            FE       = string(50);
            Appendix = '_SM_MBE_Farm_rng_15_Plim_None';
            
           
            solution_name = strcat('Solution_GA_PCL_opt_',FE,Appendix);
            sol = load(strcat(pdata.General.pathstr, filesep,'Solution',filesep, study,filesep,...
                string(pdata.General.n_wec),'_WEC',filesep, 'Irregular', filesep, 'Results', filesep,Location,...
                filesep, solution_name));

            pdata.General.n_wec     = pdata.General.n_wec;
            pdata.WEC.Radius        = sol.pdata.Results.Radius;
            pdata.WEC.RD            = sol.pdata.Results.RD;
            pdata.WEC.Draft         = sol.pdata.Results.Draft;
            pdata.WEC.PTO_damping   = sol.pdata.Results.PTO_damping(1:pdata.General.n_wec);
            pdata.WEC.PTO_stiffness = sol.pdata.Results.PTO_stiffness(1:pdata.General.n_wec);
            pdata.WEC.x             = sol.pdata.Results.x(1:pdata.General.n_wec,1);
            pdata.WEC.y             = sol.pdata.Results.y(1:pdata.General.n_wec,1);
            pdata.WEC.Distance      = sol.pdata.Results.Distance(1:pdata.General.n_wec,1);
            pdata.WEC.Angle         = sol.pdata.Results.Angle(1:pdata.General.n_wec,1);
            pdata.N_l               = 1; % number of simulations
            pdata.WaveData          = sol.pdata.WaveData;
            pdata.WaveData.Ai = sol.pdata.WaveData.Ai(1:pdata.General.n_wec,:,:);
            pdata.WaveData.DW = sol.pdata.WaveData.DW(1:pdata.General.n_wec,:,:);
            pdata.WaveData.SP = sol.pdata.WaveData.SP(1:pdata.General.n_wec,:,:);
            pdata.Opt.SatFlag       = sol.pdata.Opt.SatFlag;
            pdata.sol               = sol;

  

        catch
            error('Could not load solution for the simulation.')
        end

end


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
    rr  = pdata.WEC.Radius(1,i);


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

    pdata.WEC.Radius(idr_vec) = [];
    pdata.WEC.RD(idr_vec)     = [];
    pdata.WEC.Draft(idr_vec)  = [];

    pdata.WEC.x(:,idr_vec) = [];
    pdata.WEC.y(:,idr_vec) = [];

    pdata.WEC.Distance(:,idr_vec)  = [];
    pdata.WEC.Angle(:,idr_vec)     = [];
    pdata.WEC.PTO_damping(idr_vec,:) = [];
    pdata.WEC.PTO_stiffness(idr_vec,:) = [];

end


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function pdata = CheckDraft(pdata)

Radius = pdata.WEC.Radius;
Draft  = pdata.WEC.Draft;
RD     = pdata.WEC.RD;
x      = pdata.WEC.x;
y      = pdata.WEC.y;
Distance    = pdata.WEC.Distance;
Angle      = pdata.WEC.Angle;
K_pto = pdata.WEC.PTO_stiffness;
B_pto = pdata.WEC.PTO_damping;

id_r = [];
k = 1;
for i = 1:length(Radius)
    if Draft(i) < 0.5 || Draft(i) > 20
        id_r(k) = i;
        k = k+1;
    end
end

Radius(id_r) = [];
Draft(id_r)  = [];
RD(id_r)     = [];
x(:,id_r)      = [];
y(:,id_r)      = [];
Distance(:,id_r)   = [];
Angle(:,id_r)      = [];
K_pto(id_r,:) = [];
B_pto(id_r,:) = [];

pdata.WEC.Radius = Radius;
pdata.WEC.Draft = Draft;
pdata.WEC.RD = RD;
pdata.WEC.x = x;
pdata.WEC.y = y;
pdata.WEC.Distance = Distance;
pdata.WEC.Angle = Angle;
pdata.WEC.PTO_stiffness = K_pto;
pdata.WEC.PTO_damping = B_pto;

end

