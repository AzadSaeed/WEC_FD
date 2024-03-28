function [c, ceq] = Constr(y,pdata)

switch upper(pdata.General.CaseStudy)

    case 'P_OPT'

        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);
        np                  = pdata.Opt.np;
        xp                  = x(1:np);
        xp                  = reshape(xp,[1,2]);
        Radius              = xp(1,1);
        RD                  = xp(1,2);
        Draft               = Radius/RD;

        % Add constraints to limit draft size
        D_min               = pdata.Opt.D_min;
        D_max               = pdata.Opt.D_max;
        c_p                 = DraftConstraints(D_min, D_max, Draft);
        c                   = c_p;
        ceq                 = [];

    case 'C_OPT'

        if pdata.Opt.SatFlag
            if pdata.General.OCFlag

                error('OCFlag must be set to 0 in this release!')
                % x               = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);
                % xc              = x(1:pdata.Opt.nc);
                % c_c             = ControlSaturationConstrains(xc,pdata);
                % c               = c_c;
                % ceq             = [];

            else
                c = [];
                ceq = [];
            end
        else
            c = [];
            ceq = [];
        end

    case 'L_OPT'

        % Reverse scaling of the optimization variables
        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);

        % Separate optimization variables
        xl                  = x(1:pdata.Opt.nl);

        % Rearrange layout variables into [x, y] form
        xl                  = reshape(xl,[2,pdata.General.n_wec-1])';

        % Augment layout with the initial WEC
        xl                  = [[0,0];xl];
        Radius              = pdata.Opt.R;
        c_a                 = ArrayConstraints(xl,Radius,pdata.General.n_wec,...
            pdata.General.SafeDist);
        c                   = c_a;
        ceq                 = [];

    case 'LP_OPT'

        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);
        np                  = pdata.Opt.np;
        nl                  = pdata.Opt.nl;

        xp                  = x(1:np);
        xl                  = x((np+1):(np+nl));

        xp                  = reshape(xp,[1,2]);
        xl                  = reshape(xl,[2,pdata.General.n_wec-1])';
        xl                  = [[0,0];xl];

        Radius              = xp(1,1);
        RD                  = xp(1,2);
        Draft               = Radius/RD;

        c_a                 = ArrayConstraints(xl,Radius,pdata.General.n_wec,...
            pdata.General.SafeDist);

        D_min               = pdata.Opt.D_min;
        D_max               = pdata.Opt.D_max;
        c_p                 = DraftConstraints(D_min, D_max, Draft);

        c                   = [c_a;c_p];
        ceq                 = [];

    case 'PC_OPT'

        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);

        % Separate optimization variables
        xc                  = x(1:pdata.Opt.nc);
        xp                  = x((pdata.Opt.nc+1):(pdata.Opt.nc+pdata.Opt.np));


        % Rearrange plant variables into [Radius, Radius/Draft] form
        xp                  = reshape(xp,[1,2]);

        % Extract and assign plant variables
        Radius              = xp(1,1);
        RD                  = xp(1,2);
        Draft               = Radius/RD;

        D_min               = pdata.Opt.D_min;
        D_max               = pdata.Opt.D_max;
        c_p                 = DraftConstraints(D_min, D_max, Draft);

        if pdata.Opt.SatFlag
            if pdata.General.OCFlag

                error('OCFlag must be set to 0 in this release!')

            else
                c_c             = [];
            end
        else
            c_c             = [];
        end

        c                   = [c_p;c_c];
        ceq                 = [];

    case 'PCL_OPT'



        x                   = inv_scale_var(y,pdata.Opt.D_scale,pdata.Opt.b_scale);
        nc                  = pdata.Opt.nc;
        np                  = pdata.Opt.np;
        nl                  = pdata.Opt.nl;

        xc                  = x(1:nc);
        xp                  = x((nc+1):(nc+np));
        xl                  = x((nc+np+1):(nc+np+nl));

        xp                  = reshape(xp,[1,2]);
        xl                  = reshape(xl,[2,pdata.General.n_wec-1])';
        xl                  = [[0,0];xl];

        Radius              = xp(1,1);
        RD                  = xp(1,2);
        Draft               = Radius/RD;

        D_min               = pdata.Opt.D_min;
        D_max               = pdata.Opt.D_max;
        c_p                 = DraftConstraints(D_min, D_max, Draft);

        c_a                 = ArrayConstraints(xl,Radius,pdata.General.n_wec,...
            pdata.General.SafeDist);

        if pdata.Opt.SatFlag
            if pdata.General.OCFlag

                error('OCFlag must be set to 0 in this release!')

            else
                c_c             = [];
            end
        else
            c_c              = [];
        end

        c                   = [c_p; c_a; c_c];
        ceq                 = [];




end

end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function x = inv_scale_var(y,D,b)

y = y(:);
x = D\(y-b);

end