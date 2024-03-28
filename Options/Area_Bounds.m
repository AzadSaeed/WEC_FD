
function [lbl, ubl] = Area_Bounds(n_wec,AreaFlag)

switch AreaFlag

    case 'full'

        Location_bound               = 0.5*sqrt(20000*n_wec);
        lbl                = repmat([-Location_bound;-Location_bound],[n_wec-1,1]);
        ubl                = repmat([Location_bound;Location_bound],[n_wec-1,1]);
    case 'half'
        Location_bound               = 0.5*sqrt(20000*n_wec);
        lbl                = repmat([0;-Location_bound],[n_wec-1,1]);
        ubl                = repmat([Location_bound;Location_bound],[n_wec-1,1]);
end
end