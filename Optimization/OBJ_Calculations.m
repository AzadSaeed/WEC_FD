function varargout = OBJ_Calculations(y,pdata)

switch upper(pdata.General.CaseStudy)

    case 'P_OPT'

        % Reverse scaling of the optimization variables
        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);

        % Separate optimization variables
        xp                  = x(1:pdata.Opt.np);

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = reshape(xp,[1,2]);

        % Extract and assign plant variables
        R                       = xp(1,1);
        RD                      = xp(1,2);

        % Assign plant
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;


        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = pdata.Opt.B_pto;
        pdata.WEC.PTO_stiffness = pdata.Opt.k_pto;

        % Rearrange layout variables into [x, y] form
        pdata.WEC.x       = pdata.Opt.x;
        pdata.WEC.y       = pdata.Opt.y;



    case 'C_OPT'

        % Reverse scaling of the optimization variables
        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);

        % Separate optimization variables
        xc                  = x(1:pdata.Opt.nc);

        % Rearrange control parameters into [Damping, Stiffness]  form
        switch upper(pdata.Opt.CtrlFlag)

            % Optimize each control individually
            case 'INDIVIDUAL'

                xc                  = reshape(xc,[2,pdata.General.n_wec])';

            case 'FARM'
                xc                  = repmat(xc',[pdata.General.n_wec,1]);
        end

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = [pdata.Opt.R, pdata.Opt.RD];

        % Rearrange layout variables into [x, y] form
        xl                  = [pdata.Opt.x, pdata.Opt.y];

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = xc(:,1)';
        pdata.WEC.PTO_stiffness = xc(:,2)';

        pdata.WEC.x       = pdata.Opt.x;
        pdata.WEC.y       = pdata.Opt.y;
        pdata.WEC.Radius  = pdata.Opt.R;
        pdata.WEC.RD      = pdata.Opt.RD;
        pdata.WEC.Draft   = pdata.WEC.Radius/pdata.WEC.RD;


    case 'L_OPT'

        % Reverse scaling of the optimization variables
        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);

        % Separate optimization variables
        xl                  = x(1:pdata.Opt.nl);

        % Rearrange layout variables into [x, y] form
        xl                  = reshape(xl,[2,pdata.General.n_wec-1])';

        % Augment layout with the initial WEC
        xl                  = [[0,0];xl];

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = [pdata.Opt.R, pdata.Opt.RD];

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = pdata.Opt.B_pto;
        pdata.WEC.PTO_stiffness = pdata.Opt.k_pto;

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
        pdata.WEC.Radius   = pdata.Opt.R;
        pdata.WEC.RD       = pdata.Opt.RD;
        pdata.WEC.Draft    = pdata.WEC.Radius/pdata.WEC.RD;
        pdata.WEC.x        = xl(:,1);
        pdata.WEC.y        = xl(:,2);



    case 'LP_OPT'

        % Reverse scaling of the optimization variables
        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);

        % Separate optimization variables
        xp                  = x(1:pdata.Opt.np);
        xl                  = x((pdata.Opt.np+1):(pdata.Opt.np+pdata.Opt.nl));

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = reshape(xp,[1,2]);

        % Rearrange layout variables into [x, y] form
        xl                  = reshape(xl,[2,pdata.General.n_wec-1])';

        % Augment layout with the initial WEC
        xl                  = [[0,0];xl];

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = pdata.Opt.B_pto;
        pdata.WEC.PTO_stiffness = pdata.Opt.k_pto;

        % Extract and assign plant variables
        R                       = xp(1,1);
        RD                      = xp(1,2);

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


    case 'PC_OPT'

        % Reverse scaling of the optimization variables
        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);

        % Separate optimization variables
        xc                  = x(1:pdata.Opt.nc);
        xp                  = x((pdata.Opt.nc+1):(pdata.Opt.nc+pdata.Opt.np));

        % Rearrange control parameters into [Damping, Stiffness]  form
        switch upper(pdata.Opt.CtrlFlag)

            % Optimize each control individually
            case 'INDIVIDUAL'
                xc                  = reshape(xc,[2,pdata.General.n_wec])';

            case 'FARM'
                xc                  = repmat(xc',[pdata.General.n_wec,1]);
        end


        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = reshape(xp,[1,2]);

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = xc(:,1)';
        pdata.WEC.PTO_stiffness = xc(:,2)';

        % Extract and assign plant variables
        R                       = xp(1,1);
        RD                      = xp(1,2);


        % Rearrange layout variables into [x, y] form
        pdata.WEC.x       = pdata.Opt.x;
        pdata.WEC.y       = pdata.Opt.y;

        % Calculate HD
        pdata.WEC.Radius   = R;
        pdata.WEC.RD       = RD;
        pdata.WEC.Draft    = R/RD;



    case 'PCL_OPT'

        % Reverse scaling of the optimization variables
        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);



        % Separate optimization variables
        xc                  = x(1:pdata.Opt.nc);
        xp                  = x((pdata.Opt.nc+1):(pdata.Opt.nc+pdata.Opt.np));
        xl                  = x((pdata.Opt.nc+pdata.Opt.np+1):(pdata.Opt.nc+pdata.Opt.np+pdata.Opt.nl));


        % Rearrange layout variables into [x, y] form
        xl                  = reshape(xl,[2,pdata.General.n_wec-1])';

        % Augment layout with the initial WEC
        xl                  = [[0,0];xl];

        % Rearrange control parameters into [Damping, Stiffness]  form
        switch upper(pdata.Opt.CtrlFlag)

            % Optimize each control individually
            case 'INDIVIDUAL'
                xc                  = reshape(xc,[2,pdata.General.n_wec])';

            case 'FARM'
                xc                  = repmat(xc',[pdata.General.n_wec,1]);
        end

        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = reshape(xp,[1,2]);

        % Extract and assign PTO variables
        pdata.WEC.PTO_damping   = xc(:,1)';
        pdata.WEC.PTO_stiffness = xc(:,2)';

        % Extract and assign plant variables
        R                       = xp(1,1);
        RD                      = xp(1,2);


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

end


% P_total: annual average power from WEC farm [W]
[P_total, pdata, power_, p_mat, Pe, Pae] = WEC_Dyn(pdata);

% Power after accounting for operational availability and tranmission
% efficiency [W]
P_total_vec = P_total.*pdata.General.Op_av.*pdata.General.Trans_eff;

if strcmpi(pdata.General.CaseStudy,'SIM')
    Obj = [];
else

    switch upper(pdata.Opt.ObjForm)

        case 'POWER'

            % Power [kW]
            Power_t = sum(sum(P_total_vec))/1000;
            Obj = Power_t;

        case 'POWERPERVOLUME'

            % n_wec*Device volume [m^3]
            Volume = pdata.General.n_wec*(pi*pdata.WEC.Radius^2)*2*pdata.WEC.Draft;

            % Power [kW]
            Power_t = sum(sum(P_total_vec))/1000;

            % Energy/Volume [kW/m^3]
            Obj = -Power_t/Volume;

    end
end

if nargout == 1

    varargout  = {Obj};

elseif nargout == 2

    varargout  = {Obj, pdata};

elseif nargout >= 2

    varargout  = {Obj, pdata, P_total_vec, power_, p_mat, Pe, Pae};
end


end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function x = inv_scale_var(y,D,b)

y = y(:);
x = D\(y-b);

end
