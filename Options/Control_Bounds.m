function [lbc,ubc, nc] = Control_Bounds(n_wec,B_PTO_min,k_PTO_min,B_PTO_max,k_PTO_max, CtrlFlag)
switch upper(CtrlFlag)

    % Optimize each control individually
    case 'INDIVIDUAL'

        nc                 = 2*n_wec;
        lbc                = repmat([B_PTO_min; k_PTO_min],[n_wec,1]);
        ubc                = repmat([B_PTO_max; k_PTO_max],[n_wec,1]);

        % Optimize one control for the entire farm
    case 'FARM'
        nc                 = 2;
        lbc                = [B_PTO_min ; k_PTO_min];
        ubc                = [B_PTO_max ; k_PTO_max];
end
end