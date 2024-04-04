function [A, B, Fe] = MBE_calculations(n_wec,Layout_,Omega,Radius,RD,pdata)
% This function estimates the hydrodynamic coefficients up to second order
% using the principles of many-body expansion.
% 
% The inputs:
%            - n_wec: number of WECs in the farm
%            - Layout: farm layout where x is Layout(:,1) and y is Layout(:,2) 
%
% The outputs:
%            - A: added mass
%            - B: damping coefficient
%            - Fe: excitation force



% Original frequency corresponding to surrogate model outputs
w0     = linspace(pdata.General.wMin,pdata.General.wMax,100);

% Extract
depth  = pdata.General.depth;
x_     = Layout_(:,1);
y_     = Layout_(:,2);
xref0  = 0;

% Calculate all 1-WEC QoI
Input_1WEC = [Radius;RD];

% Create an empty SM field if not already created
if ~isfield(pdata,'SM')
    pdata.SM = [];
end


% A
a_pp0      = Calc_mean_cmt(Input_1WEC,'A',pdata.SM);
F_a_pp0    = griddedInterpolant(w0, a_pp0,'linear');
a_pp0      = F_a_pp0(Omega);

% B
b_pp0      = Calc_mean_cmt(Input_1WEC,'B',pdata.SM);
F_b_pp0    = griddedInterpolant(w0, b_pp0,'linear');
b_pp0      = F_b_pp0(Omega);


% Fe_real
fe_r_pp0   = Calc_mean_cmt(Input_1WEC,'Fe_r',pdata.SM);
F_fr_pp0   = griddedInterpolant(w0, fe_r_pp0,'linear');
fe_r_pp0   = F_fr_pp0(Omega);


% Fe_imaginary
fe_i_pp0   = Calc_mean_cmt(Input_1WEC,'Fe_i',pdata.SM);
F_fi_pp0   = griddedInterpolant(w0, fe_i_pp0,'linear');
fe_i_pp0   = F_fi_pp0(Omega);

% Fe
Fe_00      = fe_r_pp0 + 1i*fe_i_pp0;

% Calculate wave number once
wavenumber              = calcWaveNumber(Omega,depth,pdata.General.g,0);

% Number of frequencies
nfreq            = length(Omega);

if pdata.General.n_wec > 1

    % Define all coefficients in a 3D matrix, where the third dimension is
    % associated with frequency and rows and columns are associated with the
    % interaction of p and q body

    A                = zeros(n_wec,n_wec,nfreq);
    B                = zeros(n_wec,n_wec,nfreq);
    Fe               = zeros(1,n_wec,nfreq);
    Delta_a          = zeros(n_wec,n_wec,nfreq);
    Delta_b          = zeros(n_wec,n_wec,nfreq);
    Delta_fr         = zeros(n_wec,n_wec,nfreq);
    Delta_fi         = zeros(n_wec,n_wec,nfreq);
    ix_v             = 1:n_wec;

    for i = 1:n_wec

        % Define pth body
        p_body  = i;

        % Remove pth body from the set of indices
        ix_v(p_body) = [];

        % Perform this for all remaining bodies
        for j = 1:length(ix_v)

            % Define qth body
            q_body  = ix_v(j);

            % Prescribe the location of Pth body
            xc_p    = x_(p_body);
            yc_p    = y_(p_body);

            % Prescribe the location of qth body
            xc_q    = x_(q_body);
            yc_q    = y_(q_body);

            % Find their distance and angle
            delta_x  = xc_q - xc_p;
            delta_y  = yc_q - yc_p;
            l_pq     = sqrt(delta_x^2 + delta_y^2);
            theta_pq = atan2((delta_y*sign(delta_y)), delta_x);

            % Create appropriate inputs for 2-WEC surrogate models
            Input_2WEC = [Radius;RD;l_pq;theta_pq];


            % Calculate Additive effect for A_11:

            % Original frequency
            a_pp_2wec                = Calc_mean_cmt(Input_2WEC,'A_11',pdata.SM);
            F_a_pp2wec               = griddedInterpolant(w0, a_pp_2wec,'linear');
            a_pp_2wec                = F_a_pp2wec(Omega);
            delta_a                  = a_pp_2wec - a_pp0;
            Delta_a(p_body,q_body,:) = delta_a;


            % Calculate Interaction effect for A_12
            a_pq                     = Calc_mean_cmt(Input_2WEC,'A_12',pdata.SM);
            F_a_pq2wec               = griddedInterpolant(w0, a_pq,'linear');
            a_pq                     = F_a_pq2wec(Omega);
            A(p_body,q_body,:)       = a_pq;


            % Calculate Additive effect for B_11
            b_pp_2wec                = Calc_mean_cmt(Input_2WEC,'B_11',pdata.SM);
            F_b_pp2wec               = griddedInterpolant(w0, b_pp_2wec,'linear');
            b_pp_2wec                = F_b_pp2wec(Omega);
            delta_b                  = b_pp_2wec - b_pp0;
            Delta_b(p_body,q_body,:) = delta_b;

            % Calculate Interaction effect for B_12
            b_pq                     = Calc_mean_cmt(Input_2WEC,'B_12',pdata.SM);
            F_b_pq2wec               = griddedInterpolant(w0, b_pq,'linear');
            b_pq                     = F_b_pq2wec(Omega);
            B(p_body,q_body,:)       = b_pq;

            % Calculate excitation force
            fe_00_tr                = Fe_00;

            % Calculate Fe_r_11 and Fe_i_11 for pth body
            fe_r_pp_2wec             = Calc_mean_cmt(Input_2WEC,'Fe_r_11',pdata.SM);
            F_fr_pp2wec              = griddedInterpolant(w0, fe_r_pp_2wec,'linear');
            fe_r_pp_2wec             = F_fr_pp2wec(Omega);


            fe_i_pp_2wec             = Calc_mean_cmt(Input_2WEC,'Fe_i_11',pdata.SM);
            F_fi_pp2wec              = griddedInterpolant(w0, fe_i_pp_2wec,'linear');
            fe_i_pp_2wec             = F_fi_pp2wec(Omega);

            fe_2_wec                 = fe_r_pp_2wec + 1i.*fe_i_pp_2wec;

            % Calculate the distance for the additional phase shift
            Del = xc_p;

            % Calculate the additive effect
            delta_F                   = (fe_2_wec-fe_00_tr).*exp(-1i.*wavenumber*(Del));
            Delta_fr(p_body,q_body,:) = real(delta_F);
            Delta_fi(p_body,q_body,:) = imag(delta_F);

        end

        % Assign all added mass
        Delta_a_                 = sum(Delta_a(i,:,:),2);
        A(p_body,p_body,:)       = a_pp0 + squeeze(Delta_a_);

        % Assign all damping coefficients
        Delta_b_                 = sum(Delta_b(i,:,:),2);
        B(p_body,p_body,:)       = b_pp0 + squeeze(Delta_b_);

        % Assign all excitation forces
        Delta_Fr                 = sum(Delta_fr(i,:,:),2);
        Delta_Fi                 = sum(Delta_fi(i,:,:),2);
        Delta_F                  = Delta_Fr + 1i*Delta_Fi;
        Fe(1,p_body,:)           = fe_00_tr.*exp(1i.*wavenumber*(xc_p-xref0)) + squeeze(Delta_F);

        ix_v                     = 1:n_wec;

    end

else

    A       = reshape(a_pp0,[1,1,nfreq]);
    B       = reshape(b_pp0,[1,1,nfreq]);
    xc_p    = Layout_(1,1);
    Fe_     = Fe_00.*exp(1i.*wavenumber*(xc_p-xref0));
    Fe      = reshape(Fe_,[1,1,nfreq]);

end


[A_cf, B_cf, Fe_cf] = NormalizingCoeff(Radius, pdata.WEC.Draft, Omega, pdata.General.rho, pdata.General.g);
B_cf                = repmat(B_cf', n_wec*n_wec, 1);
B_cf_               = reshape(B_cf,[n_wec,n_wec,nfreq]);
A                   = A*A_cf;
B                   = B.*B_cf_;
Fe                  = Fe*Fe_cf;

end



