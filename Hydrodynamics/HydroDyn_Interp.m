function [A_eff, B_eff, Fe, A, B, Hs_k, Mass, pdata] = HydroDyn_Interp(Omega, pdata)
% This function calculates the hydrostatic and hydrodynamic coefficients
% based on the options provided by the user. If pdata.Hyrdoflag is set to
% "MS", then the code uses multi scattering and if the flag is set to
% "MBE", then the code uses many-body expansion principles along with the
% surrogate models that have been already developed.
%
% Contributor: Saeed Azad, PhD


% Estimate Mass and hydrostatic coefficients
S          = pi*(pdata.WEC.Radius)^2;
G          = pdata.General.rho*pdata.General.g*S;
Hs_k       = G;
Mass       = pdata.General.rho*S*(2*pdata.WEC.Draft);


switch upper(pdata.HydroFlag)

    case 'MS'


        pdata.WEC.w = Omega;
        [A, B, Fe]  = MS_calculations(Omega,pdata);


        % Calculate the A_eff = A(w)+M
        A_eff      = A + Mass;

        % Calculate B_eff = B(w) + PTO_damping
        B_PTO     = pdata.WEC.PTO_damping';
        B_eff     = B + B_PTO;

    case 'SM_MBE'

        Layout_    = [pdata.WEC.x,pdata.WEC.y];
        RD         = pdata.WEC.RD;

        % Estimate hydrodynamic coefficients
        [A, B, Fe] = MBE_calculations(pdata.General.n_wec, Layout_,Omega, pdata.WEC.Radius, RD, pdata);


        % Calculate the A_eff = A(w)+M
        A_eff      = A + Mass;

        % Calculate B_eff = B(w) + PTO_damping
        B_PTO     = pdata.WEC.PTO_damping';
        B_eff     = B + B_PTO;


    case 'MS_MBE'

        Layout_     = [pdata.WEC.x,pdata.WEC.y];
        RD          = pdata.WEC.RD;
        pdata.WEC.w = Omega;

        % Estimate hydrodynamic coefficients
        [A, B, Fe] = MS_MBE_calculations(pdata.General.n_wec, Layout_,Omega, pdata.WEC.Radius, RD, pdata);


        % Calculate the A_eff = A(w)+M
        A_eff      = A + Mass;

        % Calculate B_eff = B(w) + PTO_damping
        B_PTO     = pdata.WEC.PTO_damping';
        B_eff     = B + B_PTO;

end




end


