function [A, B, Fe] = MS_MBE_calculations(n_wec,Layout_,Omega,Radius,RD,pdata)
% This function estimates the hydrodynamic coefficients up to second order
% using the principles of many-body expansion, but instead of using
% surrogate models, it directly uses multi-scattering. 
% 
% The inputs:
%            - n_wec: number of WECs in the farm
%            - Layout: farm layout where x is Layout(:,1) and y is Layout(:,2) 
%
% The outputs:
%            - A: added mass
%            - B: damping coefficient
%            - Fe: excitation force


% Extract
depth  = pdata.General.depth;
x_     = Layout_(:,1);
y_     = Layout_(:,2);
xref0  = 0;

% Call Multi-scattering- for 1WEC
PDATA = pdata;
PDATA.General.n_wec = 1;
PDATA.WEC.x = 0;
PDATA.WEC.y = 0;
PDATA.WEC.Radius = Radius;
PDATA.WEC.Draft = Radius/RD;
[a_pp0_, b_pp0_, Fe_00_]  = MS_calculations(Omega,PDATA);
clear PDATA


% A
a_pp0 = squeeze(a_pp0_);

% B
b_pp0 = squeeze(b_pp0_);

% Fe
Fe_00 = squeeze(Fe_00_);

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

            % Call Multi-scattering- for 1WEC
            PDATA = pdata;
            PDATA.General.n_wec = 2;
            PDATA.WEC.Radius = Radius;
            PDATA.WEC.Draft = Radius/RD;
            PDATA.WEC.x = [0;l_pq*cos(theta_pq)];
            PDATA.WEC.y = [0;l_pq*sin(theta_pq)];
            [a_2wec, b_2wec, Fe_2wec]  = MS_calculations(Omega,PDATA);
            clear PDATA




            % Calculate Additive effect for A_11:

            % Original frequency
            a_pp_2wec                = squeeze(a_2wec(1,1,:));
            delta_a                  = a_pp_2wec - a_pp0;
            Delta_a(p_body,q_body,:) = delta_a;


            % Calculate Interaction effect for A_12
            a_pq                     = squeeze(a_2wec(1,2,:));
            A(p_body,q_body,:)       = a_pq;


            % Calculate Additive effect for B_11
            b_pp_2wec                = squeeze(b_2wec(1,1,:));
            delta_b                  = b_pp_2wec - b_pp0;
            Delta_b(p_body,q_body,:) = delta_b;

            % Calculate Interaction effect for B_12
            b_pq                     = squeeze(b_2wec(1,2,:));
            B(p_body,q_body,:)       = b_pq;

            % Calculate excitation force
            fe_00_tr                = Fe_00;

            % Calculate Fe_r_11 and Fe_i_11 for pth body
            fe_2_wec                 = squeeze(Fe_2wec(1,1,:));

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

    % needs work
    A       = reshape(a_pp0,[1,1,nfreq]);
    B       = reshape(b_pp0,[1,1,nfreq]);
    xc_p    = Layout_(1,1);
    Fe_     = Fe_00.*exp(1i.*wavenumber*(xc_p-xref0));
    Fe      = reshape(Fe_,[1,1,nfreq]);

end


end
