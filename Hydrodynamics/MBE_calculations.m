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



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

% function pdata = LoadALLSMs(pdata, flag)
% 
% switch flag
% 
%     case 'Non-CMT'
% 
%         % Load all 1-WEC models
%         ModelSet   = {'A','B','Fe_r', 'Fe_i'};
%         Model_dir  = strcat(pdata.General.pathstr, filesep,'SurrogateModels');
%         for i = 1:4
%             Name = strcat(Model_dir,filesep,'SM_1WEC_',ModelSet{1,i},'_SM_1WEC_lhs_rd_ratio');
%             load(Name)
%             pdata.SM.Net_1WEC{1,i} = net;
%         end
% 
%         % Load all 2-WEC models
%         ModelSet   = {'A_11','A_12','B_11','B_12','Fe_r_11', 'Fe_i_11'};
% 
%         for i = 1:6
%             Name = strcat(Model_dir,filesep,'SM_2WEC_',ModelSet{1,i},'_SM_2WEC_Grid_rd_ratio');
%             load(Name)
%             pdata.SM.Net_2WEC{1,i} = net;
%         end
% 
%     case 'CMT'
%         
%         Model_dir  = strcat(pdata.General.pathstr, filesep,'SurrogateModels');
% 
%         % Load all 1-WEC models
%         ModelSet   = {'A','B','Fe_r', 'Fe_i'};
% 
%         for i = 1:4
%             Name = strcat(Model_dir,filesep,'CMT_1WEC_QBC_',ModelSet{1,i});
%             load(Name)
%             pdata.SM.Net_1WEC{1,i} = Committee;
%         end
% 
%         % Load all 2-WEC models
%         ModelSet   = {'A_11','A_12','B_11','B_12','Fe_r_11', 'Fe_i_11'};
% 
%         for i = 1:6
%             Name = strcat(Model_dir,filesep,'CMT_2WEC_QBC_',ModelSet{1,i});
%             load(Name)
%             pdata.SM.Net_2WEC{1,i} = Committee;
%         end
% 
% end
% 
% end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% 
% function Create2figs(A,B,Fe,pdata)
% 
% Omega = pdata.w;
% 
% ValidateMBE('2_wec_A_11',squeeze(A(1,1,:)),pdata.MS.a_11,Omega,pdata.plot_dir,'$A_{11}$')
% ValidateMBE('2_wec_A_12',squeeze(A(1,2,:)),pdata.MS.a_12,Omega,pdata.plot_dir,'$A_{12}$')
% ValidateMBE('2_wec_A_21',squeeze(A(2,1,:)),pdata.MS.a_21,Omega,pdata.plot_dir,'$A_{21}$')
% ValidateMBE('2_wec_A_22',squeeze(A(2,2,:)),pdata.MS.a_22,Omega,pdata.plot_dir,'$A_{22}$')
% 
% ValidateMBE('2_wec_B_11',squeeze(B(1,1,:)),pdata.MS.b_11,Omega,pdata.plot_dir,'$B_{11}$')
% ValidateMBE('2_wec_B_12',squeeze(B(1,2,:)),pdata.MS.b_12,Omega,pdata.plot_dir,'$B_{12}$')
% ValidateMBE('2_wec_B_21',squeeze(B(2,1,:)),pdata.MS.b_21,Omega,pdata.plot_dir,'$B_{21}$')
% ValidateMBE('2_wec_B_22',squeeze(B(2,2,:)),pdata.MS.b_22,Omega,pdata.plot_dir,'$B_{22}$')
% 
% ValidateMBE('2_wec_Fe_r_11',squeeze(real(Fe(1,1,:))),pdata.MS.fe_r_11,Omega,pdata.plot_dir,'$Fe_{r_{11}}$')
% ValidateMBE('2_wec_Fe_i_11',squeeze(imag(Fe(1,1,:))),pdata.MS.fe_i_11,Omega,pdata.plot_dir,'$Fe_{i_{11}}$')
% ValidateMBE('2_wec_Fe_r_22',squeeze(real(Fe(1,2,:))),pdata.MS.fe_r_22,Omega,pdata.plot_dir,'$Fe_{r_{22}}$')
% ValidateMBE('2_wec_Fe_i_22',squeeze(imag(Fe(1,2,:))),pdata.MS.fe_i_22,Omega,pdata.plot_dir,'$Fe_{i_{22}}$')
% 
% end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

% function Create3figs(A,B,Fe,pdata)
% 
% Omega = pdata.w;
% 
% ValidateMBE('3_wec_A_11',squeeze(A(1,1,:)),pdata.MS.a_11,Omega,pdata.plot_dir,'$A_{11}$')
% ValidateMBE('3_wec_A_12',squeeze(A(1,2,:)),pdata.MS.a_12,Omega,pdata.plot_dir,'$A_{12}$')
% ValidateMBE('3_wec_A_13',squeeze(A(1,3,:)),pdata.MS.a_13,Omega,pdata.plot_dir,'$A_{13}$')
% ValidateMBE('3_wec_A_21',squeeze(A(2,1,:)),pdata.MS.a_21,Omega,pdata.plot_dir,'$A_{21}$')
% ValidateMBE('3_wec_A_22',squeeze(A(2,2,:)),pdata.MS.a_22,Omega,pdata.plot_dir,'$A_{22}$')
% ValidateMBE('3_wec_A_23',squeeze(A(2,3,:)),pdata.MS.a_23,Omega,pdata.plot_dir,'$A_{23}$')
% ValidateMBE('3_wec_A_31',squeeze(A(3,1,:)),pdata.MS.a_31,Omega,pdata.plot_dir,'$A_{31}$')
% ValidateMBE('3_wec_A_32',squeeze(A(3,2,:)),pdata.MS.a_32,Omega,pdata.plot_dir,'$A_{32}$')
% ValidateMBE('3_wec_A_33',squeeze(A(3,3,:)),pdata.MS.a_33,Omega,pdata.plot_dir,'$A_{33}$')
% 
% ValidateMBE('3_wec_B_11',squeeze(B(1,1,:)),pdata.MS.b_11,Omega,pdata.plot_dir,'$B_{11}$')
% ValidateMBE('3_wec_B_12',squeeze(B(1,2,:)),pdata.MS.b_12,Omega,pdata.plot_dir,'$B_{12}$')
% ValidateMBE('3_wec_B_13',squeeze(B(1,3,:)),pdata.MS.b_13,Omega,pdata.plot_dir,'$B_{13}$')
% ValidateMBE('3_wec_B_21',squeeze(B(2,1,:)),pdata.MS.b_21,Omega,pdata.plot_dir,'$B_{21}$')
% ValidateMBE('3_wec_B_22',squeeze(B(2,2,:)),pdata.MS.b_22,Omega,pdata.plot_dir,'$B_{22}$')
% ValidateMBE('3_wec_B_23',squeeze(B(2,3,:)),pdata.MS.b_23,Omega,pdata.plot_dir,'$B_{23}$')
% ValidateMBE('3_wec_B_31',squeeze(B(3,1,:)),pdata.MS.b_31,Omega,pdata.plot_dir,'$B_{31}$')
% ValidateMBE('3_wec_B_32',squeeze(B(3,2,:)),pdata.MS.b_32,Omega,pdata.plot_dir,'$B_{32}$')
% ValidateMBE('3_wec_B_33',squeeze(B(3,3,:)),pdata.MS.b_33,Omega,pdata.plot_dir,'$B_{33}$')
% 
% ValidateMBE('3_wec_Fe_r_11',squeeze(real(Fe(1,1,:))),pdata.MS.fe_r_11,Omega,pdata.plot_dir,'$Fe_{r_{11}}$')
% ValidateMBE('3_wec_Fe_i_11',squeeze(imag(Fe(1,1,:))),pdata.MS.fe_i_11,Omega,pdata.plot_dir,'$Fe_{i_{11}}$')
% ValidateMBE('3_wec_Fe_r_22',squeeze(real(Fe(1,2,:))),pdata.MS.fe_r_22,Omega,pdata.plot_dir,'$Fe_{r_{22}}$')
% ValidateMBE('3_wec_Fe_i_22',squeeze(imag(Fe(1,2,:))),pdata.MS.fe_i_22,Omega,pdata.plot_dir,'$Fe_{i_{22}}$')
% ValidateMBE('3_wec_Fe_r_33',squeeze(real(Fe(1,3,:))),pdata.MS.fe_r_33,Omega,pdata.plot_dir,'$Fe_{r_{33}}$')
% ValidateMBE('3_wec_Fe_i_33',squeeze(imag(Fe(1,3,:))),pdata.MS.fe_i_33,Omega,pdata.plot_dir,'$Fe_{i_{33}}$')
% 
% end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

% function Create4figs(A,B,Fe,pdata)
% 
% Omega = pdata.w;
% 
% ValidateMBE('4_wec_A_11',squeeze(A(1,1,:)),pdata.MS.a_11,Omega,pdata.plot_dir,'$A_{11}$')
% ValidateMBE('4_wec_A_12',squeeze(A(1,2,:)),pdata.MS.a_12,Omega,pdata.plot_dir,'$A_{12}$')
% ValidateMBE('4_wec_A_13',squeeze(A(1,3,:)),pdata.MS.a_13,Omega,pdata.plot_dir,'$A_{13}$')
% ValidateMBE('4_wec_A_14',squeeze(A(1,4,:)),pdata.MS.a_14,Omega,pdata.plot_dir,'$A_{14}$')
% ValidateMBE('4_wec_A_21',squeeze(A(2,1,:)),pdata.MS.a_21,Omega,pdata.plot_dir,'$A_{21}$')
% ValidateMBE('4_wec_A_22',squeeze(A(2,2,:)),pdata.MS.a_22,Omega,pdata.plot_dir,'$A_{22}$')
% ValidateMBE('4_wec_A_23',squeeze(A(2,3,:)),pdata.MS.a_23,Omega,pdata.plot_dir,'$A_{23}$')
% ValidateMBE('4_wec_A_24',squeeze(A(2,4,:)),pdata.MS.a_24,Omega,pdata.plot_dir,'$A_{24}$')
% ValidateMBE('4_wec_A_31',squeeze(A(3,1,:)),pdata.MS.a_31,Omega,pdata.plot_dir,'$A_{31}$')
% ValidateMBE('4_wec_A_32',squeeze(A(3,2,:)),pdata.MS.a_32,Omega,pdata.plot_dir,'$A_{32}$')
% ValidateMBE('4_wec_A_33',squeeze(A(3,3,:)),pdata.MS.a_33,Omega,pdata.plot_dir,'$A_{33}$')
% ValidateMBE('4_wec_A_34',squeeze(A(3,4,:)),pdata.MS.a_34,Omega,pdata.plot_dir,'$A_{34}$')
% ValidateMBE('4_wec_A_41',squeeze(A(4,1,:)),pdata.MS.a_41,Omega,pdata.plot_dir,'$A_{41}$')
% ValidateMBE('4_wec_A_42',squeeze(A(4,2,:)),pdata.MS.a_42,Omega,pdata.plot_dir,'$A_{42}$')
% ValidateMBE('4_wec_A_43',squeeze(A(4,3,:)),pdata.MS.a_43,Omega,pdata.plot_dir,'$A_{43}$')
% ValidateMBE('4_wec_A_44',squeeze(A(4,4,:)),pdata.MS.a_44,Omega,pdata.plot_dir,'$A_{44}$')
% 
% 
% ValidateMBE('4_wec_B_11',squeeze(B(1,1,:)),pdata.MS.b_11,Omega,pdata.plot_dir,'$B_{11}$')
% ValidateMBE('4_wec_B_12',squeeze(B(1,2,:)),pdata.MS.b_12,Omega,pdata.plot_dir,'$B_{12}$')
% ValidateMBE('4_wec_B_13',squeeze(B(1,3,:)),pdata.MS.b_13,Omega,pdata.plot_dir,'$B_{13}$')
% ValidateMBE('4_wec_B_14',squeeze(B(1,4,:)),pdata.MS.b_14,Omega,pdata.plot_dir,'$B_{14}$')
% ValidateMBE('4_wec_B_21',squeeze(B(2,1,:)),pdata.MS.b_21,Omega,pdata.plot_dir,'$B_{21}$')
% ValidateMBE('4_wec_B_22',squeeze(B(2,2,:)),pdata.MS.b_22,Omega,pdata.plot_dir,'$B_{22}$')
% ValidateMBE('4_wec_B_23',squeeze(B(2,3,:)),pdata.MS.b_23,Omega,pdata.plot_dir,'$B_{23}$')
% ValidateMBE('4_wec_B_24',squeeze(B(2,4,:)),pdata.MS.b_24,Omega,pdata.plot_dir,'$B_{24}$')
% ValidateMBE('4_wec_B_31',squeeze(B(3,1,:)),pdata.MS.b_31,Omega,pdata.plot_dir,'$B_{31}$')
% ValidateMBE('4_wec_B_32',squeeze(B(3,2,:)),pdata.MS.b_32,Omega,pdata.plot_dir,'$B_{32}$')
% ValidateMBE('4_wec_B_33',squeeze(B(3,3,:)),pdata.MS.b_33,Omega,pdata.plot_dir,'$B_{33}$')
% ValidateMBE('4_wec_B_34',squeeze(B(3,4,:)),pdata.MS.b_34,Omega,pdata.plot_dir,'$B_{34}$')
% ValidateMBE('4_wec_B_41',squeeze(B(4,1,:)),pdata.MS.b_41,Omega,pdata.plot_dir,'$B_{41}$')
% ValidateMBE('4_wec_B_42',squeeze(B(4,2,:)),pdata.MS.b_42,Omega,pdata.plot_dir,'$B_{42}$')
% ValidateMBE('4_wec_B_43',squeeze(B(4,3,:)),pdata.MS.b_43,Omega,pdata.plot_dir,'$B_{43}$')
% ValidateMBE('4_wec_B_44',squeeze(B(4,4,:)),pdata.MS.b_44,Omega,pdata.plot_dir,'$B_{44}$')
% 
% ValidateMBE('4_wec_Fe_r_11',squeeze(real(Fe(1,1,:))),pdata.MS.fe_r_11,Omega,pdata.plot_dir,'$Fe_{r_{11}}$')
% ValidateMBE('4_wec_Fe_i_11',squeeze(imag(Fe(1,1,:))),pdata.MS.fe_i_11,Omega,pdata.plot_dir,'$Fe_{i_{11}}$')
% ValidateMBE('4_wec_Fe_r_22',squeeze(real(Fe(1,2,:))),pdata.MS.fe_r_22,Omega,pdata.plot_dir,'$Fe_{r_{22}}$')
% ValidateMBE('4_wec_Fe_i_22',squeeze(imag(Fe(1,2,:))),pdata.MS.fe_i_22,Omega,pdata.plot_dir,'$Fe_{i_{22}}$')
% ValidateMBE('4_wec_Fe_r_33',squeeze(real(Fe(1,3,:))),pdata.MS.fe_r_33,Omega,pdata.plot_dir,'$Fe_{r_{33}}$')
% ValidateMBE('4_wec_Fe_i_33',squeeze(imag(Fe(1,3,:))),pdata.MS.fe_i_33,Omega,pdata.plot_dir,'$Fe_{i_{33}}$')
% ValidateMBE('4_wec_Fe_r_44',squeeze(real(Fe(1,4,:))),pdata.MS.fe_r_44,Omega,pdata.plot_dir,'$Fe_{r_{44}}$')
% ValidateMBE('4_wec_Fe_i_44',squeeze(imag(Fe(1,4,:))),pdata.MS.fe_i_44,Omega,pdata.plot_dir,'$Fe_{i_{44}}$')
% 
% end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


% function ValidateMBE(flag,Y,MS,Omega,plot_dir,Ylab)
% 
% 
% % Interpolate MBE solution at 200 points
% xq    = linspace(0.1,7,200);
% Y     = interp1(Omega,Y,xq);
% 
% 
% commonFigureSetup;
% commonFigureProperties;
% hf= figure;
% plot(xq,Y,'-','linewidth',linewidth,'markersize',14,'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(9,:));
% hold on
% plot(Omega,MS,'.','linewidth',linewidth,'markersize',14,'Color',C.red(9,:),'MarkerEdgeColor',C.red(9,:));
% xlabel('$\omega$'); ylabel(Ylab);
% savename = strcat(plot_dir,filesep,'MBE_Validation_',flag,'.pdf');
% exportgraphics(hf,savename)
% close(hf)
% 
% end
