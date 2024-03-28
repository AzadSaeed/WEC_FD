function pdata = PostProcess(y,fval,pdata)


% Reverse scaling of the optimization variables
x = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);



switch upper(pdata.General.CaseStudy)

    case 'P_OPT'

        % Extract optimization variables
        xp                     = x(1:pdata.Opt.np);

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = pdata.Opt.B_pto;
        pdata.WEC.PTO_stiffness = pdata.Opt.k_pto;

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp  = reshape(xp,[1,2]);
        R   = xp(1,1);
        RD  = xp(1,2);

        % Rearrange layout variables into [x, y] form
        xl                  = [pdata.Opt.x, pdata.Opt.y];


        % Define the reference WEC
        ref                     = xl(1,:);

        % Find locations wrt [0, 0]
        Dist(1,1)           = 0;
        Angle(1,1)          = 0;

        for i=2:pdata.General.n_wec
            xl1             = xl;
            delta_x         = xl1(i,1) - ref(1,1);
            delta_y         = xl1(i,2) - ref(1,2);
            Dist(i,1)       = sqrt(delta_x.^2 + delta_y.^2);
            Angle(i,1)      = atan2(delta_y,delta_x);
        end

        pdata.WEC.Distance = Dist;
        pdata.WEC.Angle    = Angle;
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;
        pdata.WEC.x        = xl(:,1);
        pdata.WEC.y        = xl(:,2);

        pdata.Results.Distance        = Dist;
        pdata.Results.Angle           = Angle;
        pdata.Results.Radius          = R;
        pdata.Results.RD              = RD;
        pdata.Results.Draft           = R/RD;
        pdata.Results.x               = xl(:,1);
        pdata.Results.y               = xl(:,2);
        pdata.Results.PTO_damping     = pdata.Opt.B_pto;
        pdata.Results.PTO_stiffness   = pdata.Opt.k_pto;
        pdata.Results.Objective       = fval;


        pdata = sim_opt(pdata);



    case 'C_OPT'


        % Extract optimization variables
        xc                  = x(1:pdata.Opt.nc);

        % Rearrange control parameters into [Damping, Stiffness]  form
        switch upper(pdata.Opt.CtrlFlag)
            % Optimize each control individually
            case 'INDIVIDUAL'
                xc                  = reshape(xc,[2,pdata.General.n_wec])';
                % Optimize one set of control parameters for the farm
            case 'FARM'
                xc                  = repmat(xc',[pdata.General.n_wec,1]);
        end

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = [pdata.Opt.R, pdata.Opt.RD];
        R                       = xp(1,1);
        RD                      = xp(1,2);


        % Rearrange layout variables into [x, y] form
        xl                  = [pdata.Opt.x, pdata.Opt.y];

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = xc(:,1)';
        pdata.WEC.PTO_stiffness = xc(:,2)';

        % Define the reference WEC
        ref                     = xl(1,:);

        % Find locations wrt [0, 0]
        Dist(1,1)           = 0;
        Angle(1,1)          = 0;

        for i=2:pdata.General.n_wec

            xl1             = xl;
            delta_x         = xl1(i,1) - ref(1,1);
            delta_y         = xl1(i,2) - ref(1,2);
            Dist(i,1)       = sqrt(delta_x.^2 + delta_y.^2);
            Angle(i,1)      = atan2(delta_y,delta_x);

        end

        % Calculate HD
        pdata.WEC.Distance = Dist;
        pdata.WEC.Angle    = Angle;
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;
        pdata.WEC.x        = xl(:,1);
        pdata.WEC.y        = xl(:,2);

        pdata.Results.Distance        = Dist;
        pdata.Results.Angle           = Angle;
        pdata.Results.Radius          = R;
        pdata.Results.RD              = RD;
        pdata.Results.Draft           = R/RD;
        pdata.Results.x               = xl(:,1);
        pdata.Results.y               = xl(:,2);
        pdata.Results.PTO_damping     = xc(:,1)';
        pdata.Results.PTO_stiffness   = xc(:,2)';
        pdata.Results.Objective       = fval;

        pdata = sim_opt(pdata);



    case 'L_OPT'


        % Extract optimization variables
        xl                      = x(1:pdata.Opt.nl);

        % Rearrange layout variables into [x, y] form
        xl  = reshape(xl,[2,pdata.General.n_wec-1])';

        % Augment layout with the initial WEC
        xl  = [[0,0];xl];

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = pdata.Opt.B_pto;
        pdata.WEC.PTO_stiffness = pdata.Opt.k_pto;


        % Extract and assign plant variables
        R                       = pdata.Opt.R;
        RD                      = pdata.Opt.RD;

        % Define the reference WEC
        ref                     = xl(1,:);

        % Find locations wrt [0, 0]
        Dist(1,1)           = 0;
        Angle(1,1)          = 0;

        for i=2:pdata.General.n_wec
            xl1             = xl;
            delta_x         = xl1(i,1) - ref(1,1);
            delta_y         = xl1(i,2) - ref(1,2);
            Dist(i,1)       = sqrt(delta_x.^2 + delta_y.^2);
            Angle(i,1)      = atan2(delta_y,delta_x);
        end

        pdata.WEC.Distance = Dist;
        pdata.WEC.Angle    = Angle;
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;
        pdata.WEC.x        = xl(:,1);
        pdata.WEC.y        = xl(:,2);

        pdata.Results.Distance        = Dist;
        pdata.Results.Angle           = Angle;
        pdata.Results.Radius          = R;
        pdata.Results.RD              = RD;
        pdata.Results.Draft           = R/RD;
        pdata.Results.x               = xl(:,1);
        pdata.Results.y               = xl(:,2);
        pdata.Results.PTO_damping     = pdata.Opt.B_pto;
        pdata.Results.PTO_stiffness   = pdata.Opt.k_pto;
        pdata.Results.Objective       = fval;

        pdata = sim_opt(pdata);



    case 'LP_OPT'

        np                  = pdata.Opt.np;
        nl                  = pdata.Opt.nl;

        xp                  = x(1:np);
        xl                  = x((np+1):(np+nl));

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp  = reshape(xp,[1,2]);
        R   = xp(1,1);
        RD  = xp(1,2);

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = pdata.Opt.B_pto;
        pdata.WEC.PTO_stiffness = pdata.Opt.k_pto;

        % Rearrange layout variables into [x, y] form
        xl  = reshape(xl,[2,pdata.General.n_wec-1])';

        % Augment layout with the initial WEC
        xl  = [[0,0];xl];

        % Define the reference WEC
        ref                     = xl(1,:);

        % Find locations wrt [0, 0]
        Dist(1,1)           = 0;
        Angle(1,1)          = 0;

        for i=2:pdata.General.n_wec
            xl1             = xl;
            delta_x         = xl1(i,1) - ref(1,1);
            delta_y         = xl1(i,2) - ref(1,2);
            Dist(i,1)       = sqrt(delta_x.^2 + delta_y.^2);
            Angle(i,1)      = atan2(delta_y,delta_x);
        end

        pdata.WEC.Distance = Dist;
        pdata.WEC.Angle    = Angle;
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;
        pdata.WEC.x        = xl(:,1);
        pdata.WEC.y        = xl(:,2);

        pdata.Results.Distance        = Dist;
        pdata.Results.Angle           = Angle;
        pdata.Results.Radius          = R;
        pdata.Results.RD              = RD;
        pdata.Results.Draft           = R/RD;
        pdata.Results.x               = xl(:,1);
        pdata.Results.y               = xl(:,2);
        pdata.Results.PTO_damping     = pdata.Opt.B_pto;
        pdata.Results.PTO_stiffness   = pdata.Opt.k_pto;
        pdata.Results.Objective       = fval;

        pdata = sim_opt(pdata);



    case 'PC_OPT'

        % Extract optimization variables
        xc                     = x(1:pdata.Opt.nc);
        xp                     = x((pdata.Opt.nc+1):(pdata.Opt.nc+pdata.Opt.np));

        % Rearrange control parameters into [Damping, Stiffness]  form
        switch upper(pdata.Opt.CtrlFlag)

            % Optimize each control individually
            case 'INDIVIDUAL'
                xc                  = reshape(xc,[2,pdata.General.n_wec])';
                % Optimize one set of control parameters for the farm
            case 'FARM'
                xc                  = repmat(xc',[pdata.General.n_wec,1]);
        end

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp  = reshape(xp,[1,2]);
        R   = xp(1,1);
        RD  = xp(1,2);

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = xc(:,1)';
        pdata.WEC.PTO_stiffness = xc(:,2)';

        % Rearrange layout variables into [x, y] form
        xl                  = [pdata.Opt.x, pdata.Opt.y];

        % Define the reference WEC
        ref                     = xl(1,:);

        % Find locations wrt [0, 0]
        Dist(1,1)           = 0;
        Angle(1,1)          = 0;

        for i=2:pdata.General.n_wec
            xl1             = xl;
            delta_x         = xl1(i,1) - ref(1,1);
            delta_y         = xl1(i,2) - ref(1,2);
            Dist(i,1)       = sqrt(delta_x.^2 + delta_y.^2);
            Angle(i,1)      = atan2(delta_y,delta_x);
        end

        pdata.WEC.Distance = Dist;
        pdata.WEC.Angle    = Angle;
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;
        pdata.WEC.x        = xl(:,1);
        pdata.WEC.y        = xl(:,2);


        pdata.Results.Distance        = Dist;
        pdata.Results.Angle           = Angle;
        pdata.Results.Radius          = R;
        pdata.Results.RD              = RD;
        pdata.Results.Draft           = R/RD;
        pdata.Results.x               = xl(:,1);
        pdata.Results.y               = xl(:,2);
        pdata.Results.PTO_damping     = pdata.WEC.PTO_damping;
        pdata.Results.PTO_stiffness   = pdata.WEC.PTO_stiffness;
        pdata.Results.Objective       = fval;

        pdata = sim_opt(pdata);


    case 'PCL_OPT'




        % Extract optimization variables
        xc                     = x(1:pdata.Opt.nc);
        xp                     = x((pdata.Opt.nc+1):(pdata.Opt.nc+pdata.Opt.np));
        xl                     = x((pdata.Opt.nc+pdata.Opt.np+1):(pdata.Opt.nc+pdata.Opt.np+pdata.Opt.nl));

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp  = reshape(xp,[1,2]);
        R   = xp(1,1);
        RD  = xp(1,2);


        % Rearrange layout variables into [x, y] form
        xl  = reshape(xl,[2,pdata.General.n_wec-1])';

        % Augment layout with the initial WEC
        xl  = [[0,0];xl];

        % Rearrange control parameters into [Damping, Stiffness]  form
        switch upper(pdata.Opt.CtrlFlag)

            % Optimize each control individually
            case 'INDIVIDUAL'
                xc                  = reshape(xc,[2,pdata.General.n_wec])';

                % Optimize one set of control parameters for the farm
            case 'FARM'
                xc                  = repmat(xc',[pdata.General.n_wec,1]);

        end

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = xc(:,1)';
        pdata.WEC.PTO_stiffness = xc(:,2)';




        % Define the reference WEC
        ref = xl(1,:);

        % Find locations wrt [0, 0]
        Dist(1,1)           = 0;
        Angle(1,1)          = 0;

        for i=2:pdata.General.n_wec

            xl1             = xl;
            delta_x         = xl1(i,1) - ref(1,1);
            delta_y         = xl1(i,2) - ref(1,2);
            Dist(i,1)       = sqrt(delta_x.^2 + delta_y.^2);
            Angle(i,1)      = atan2(delta_y,delta_x);

        end

        pdata.WEC.Distance = Dist;
        pdata.WEC.Angle    = Angle;
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;
        pdata.WEC.x        = xl(:,1);
        pdata.WEC.y        = xl(:,2);

        pdata.Results.Distance        = Dist;
        pdata.Results.Angle           = Angle;
        pdata.Results.Radius          = R;
        pdata.Results.RD              = RD;
        pdata.Results.Draft           = R/RD;
        pdata.Results.x               = xl(:,1);
        pdata.Results.y               = xl(:,2);
        pdata.Results.PTO_damping     = xc(:,1)';
        pdata.Results.PTO_stiffness   = xc(:,2)';
        pdata.Results.Objective       = fval;

        pdata = sim_opt(pdata);



end



% Remove surrogate models to reduce size of solution
if isfield(pdata,'SM')
    pdata = rmfield(pdata,'SM');
end


% Create directory for locations
if strcmpi(pdata.WaveData.Type,'Irregular')
    Lo_dir = strcat(pdata.General.Solution_dir,filesep,pdata.WaveData.Region);
else
    Lo_dir = strcat(pdata.General.Solution_dir);

end



if ~isfolder(Lo_dir)
    mkdir(Lo_dir)
end



% Scaled, randomized starting point
if ~isfield(pdata.Opt,'rng')
    ss = '5';
else
    ss = string(pdata.Opt.rng);
end


% Set up naming for power limits
if pdata.Opt.SatFlag ==0

    Pname = '_Plim_None';

elseif pdata.Opt.SatFlag ==1

    Pname = strcat('_Plim_',string(pdata.General.PowerLimit));

end


% Optimization-specific name
switch upper(pdata.Opt.Solver)

    case 'FMINCON'
        opt_ex = string(pdata.Opt.MaxFuncEval);
    case 'GA'
        opt_ex = string(pdata.Opt.MG);
    case 'SURROGATEOPT'
        opt_ex = string(pdata.Opt.MaxFuncEval);
end

% Save Output
if ~pdata.Opt.Hybrid

    switch upper(pdata.General.CaseStudy)

        case 'P_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',pdata.Opt.Solver,...
                '_',pdata.General.CaseStudy,'_',opt_ex,...
                '_', pdata.HydroFlag,'_rng_',ss,'_',lower(pdata.General.Location_Flag),Pname),'pdata');

        case 'C_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',pdata.Opt.Solver,...
                '_',pdata.General.CaseStudy,'_',opt_ex,...
                '_', pdata.HydroFlag,'_',pdata.Opt.CtrlFlag,'_rng_',ss,'_',...
                lower(pdata.General.Location_Flag),Pname),'pdata');

        case 'L_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',pdata.Opt.Solver,'_',...
                pdata.General.CaseStudy,'_',opt_ex,'_',...
                pdata.HydroFlag,'_rng_',ss,Pname),'pdata');

        case 'LP_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',pdata.Opt.Solver,'_',...
                pdata.General.CaseStudy,'_',opt_ex,'_',...
                pdata.HydroFlag,'_rng_',ss,Pname),'pdata');

        case 'PC_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',pdata.Opt.Solver,'_',...
                pdata.General.CaseStudy,'_',opt_ex,'_',...
                pdata.HydroFlag,'_',pdata.Opt.CtrlFlag,'_rng_',ss,'_',...
                lower(pdata.General.Location_Flag),Pname),'pdata');

        case 'PCL_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',...
                pdata.Opt.Solver,'_',pdata.General.CaseStudy,'_',...
                opt_ex,'_', pdata.HydroFlag,'_',...
                pdata.Opt.CtrlFlag,'_rng_',ss,Pname),'pdata');

    end

elseif pdata.Opt.Hybrid

    switch upper(pdata.General.CaseStudy)

        case 'L_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',pdata.Opt.Solver,'_',...
                pdata.General.CaseStudy,'_',opt_ex,'_',...
                pdata.HydroFlag,'_Hybrid',Pname),'pdata');

        case 'LP_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',pdata.Opt.Solver,'_',...
                pdata.General.CaseStudy,'_',opt_ex,'_',...
                pdata.HydroFlag,'_Hybrid',Pname),'pdata');

        case 'PCL_OPT'

            save(strcat(Lo_dir,filesep,'Solution_',...
                pdata.Opt.Solver,'_',pdata.General.CaseStudy,'_',...
                opt_ex,'_', pdata.HydroFlag,'_',pdata.Opt.CtrlFlag,...
                '_Hybrid',Pname),'pdata');

    end
end

end




%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function x = inv_scale_var(y,D,b)

y = y(:);
x = D\(y-b);

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = sim_opt(pdata)

Main_study = pdata.General.CaseStudy;
HydroMain  = pdata.HydroFlag;
pdata.General.CaseStudy = 'Sim';



% Run the simulation and calculate q-factor with multi scattering
pdata.HydroFlag = 'MS';

% Perform simulation
pdata                         = run_simulation(pdata);

% Calculate q-factor
[Q_factor, q_factor]          = calculate_q_factor(pdata);
pdata.Results.Q_factor_Ms     = Q_factor;
pdata.Results.q_factor_Ms     = q_factor;

pdata.Results.T.q_factor = q_factor;

% Print results table on the screen
pdata.Results.T


pdata.General.CaseStudy = Main_study;
pdata.HydroFlag         = HydroMain;
end



