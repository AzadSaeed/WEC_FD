function [power_sum, pdata] = FD_Model(pdata)
% This functions implements a frequency-domain analysis of the WEC 
% devices based on the prescribed dimensions, PTO parameters, and the
% user-defined hydrodynamic calculation flag in regular or irregular 
% waves.
% 
%-------------------------------------------------------------------------%
% Forces = F_Hydrodynamic + F_Hydrostatic +  F_PTO
%        = F_Excitation + F_Radiation + F_Hydrostatic +  F_PTO
%        = F_FK + F_scattering(diffraction) + F_Radiation + F_Hydrostatic 
%          +  F_PTO
%
%         Excitation force  = Froude-Krylov + Diffraction (or
%                             scattering force)
%
%         Radiation force   = Radiation damping + Added mass
%                           = - (-w^2*A(w) + j*w*B(w))Z(w)
%                              where Z(w) is the heave motion, A is
%                              added mass, B is radiation damping
%
%         Hydrostatic force = Buoyancy force
%                           = (F_hs(w) = - KZ(w))
%
%         PTO force         = -PTO.damping*dZ/dt - PTO.stiffness*Z
%
%         Units used in Nemoh (from Nemoh's documentation)
%         Added mass unit          = kg
%         Force unit               = N
%         Moment unit              = Nm
%         Damping coefficient unit = kg/s
%         Normalized force unit    = N/m
%         Normalized moment        = N
%         RAO_translation          = m/m
%         RAO_rotation             = deg/m
%
% Contributor: Saeed Azad, PhD
%-------------------------------------------------------------------------% 

n_wec = pdata.General.n_wec;

switch upper(pdata.WaveData.Type)

    case 'REGULAR'

        
        Hs = pdata.WaveData.Hs_regular;
        Te = pdata.WaveData.Te_regular;
        Omega = (2*pi)/Te;
        nw = length(Omega);
        [A_eff,~,Fe,Hd_A,Hd_B,Hs_K,Mass,pdata] = HydroDyn_Interp(Omega,pdata);

        % Rearrange and squeeze Fe
        if pdata.General.n_wec > 1
            Fe_ = squeeze(reshape(Fe,[n_wec,1,nw]));
        else
            Fe_ = squeeze(reshape(Fe,[n_wec,1,nw]))';
        end

        Fex_ = pdata.WaveData.S.amplitude*Fe_;
        pdata.WaveData.nodes = [Hs,Te];
        pdata.WaveData.nr = 1;


        if pdata.General.FindOptimalPTO
            pdata.WEC.PTO.k_optimal     = (Omega^2.*(sum(A_eff,2))- Hs_K)';
            pdata.WEC.PTO.B_optimal     = (sum(Hd_B,2))';
            pdata.Opt.B_pto             = pdata.WEC.PTO.B_optimal;
            pdata.Opt.k_pto             = pdata.WEC.PTO.k_optimal;
            pdata.WEC.PTO_stiffness     = pdata.WEC.PTO.k_optimal;
            pdata.WEC.PTO_damping       = pdata.WEC.PTO.B_optimal;
        end



    case 'IRREGULAR'

        % Discretize angular frequency at nr points (for the superposition
        % of reglar waves)
        Omega = (pdata.General.wMin:(pdata.General.wMax-pdata.General.wMin...
            )/(pdata.General.nbfreq-1):pdata.General.wMax)';

        % Estimate the hydrodynamic coefficients
        [~,~,Fe,Hd_A,Hd_B,Hs_K,Mass,pdata] = HydroDyn_Interp(Omega,pdata);

        % Rearrange and squeeze Fe
        nw = length(Omega);
        if pdata.General.n_wec > 1
            Fe_ = squeeze(reshape(Fe,[n_wec,1,nw]));
        else
            Fe_ = squeeze(reshape(Fe,[n_wec,1,nw]))';
        end

        % preallocate for faster computations
        nr = pdata.WaveData.nr;
        Fex_ = NaN(n_wec,nr, length(pdata.WaveData.nodes(:,1)));
        
        % Multiply the excitation force by wave amplitude
        for ii = 1:length(pdata.WaveData.nodes(:,1))
            Fex_(:,:,ii) = pdata.WaveData.Ai(:,:,ii).*Fe_;
        end

end

% Prescribe PTO parameters
PTO_damping = pdata.WEC.PTO_damping';
PTO_stiffness = pdata.WEC.PTO_stiffness';

% Manipulate matrix dimensions based on the number of devices
if pdata.General.n_wec > 1
    % Sum elements of Added mass and damping coefficient
    Hd_B_            = squeeze(sum(Hd_B,2));
    Hd_A_            = squeeze(sum(Hd_A,2));
else
    % Sum elements of Added mass and damping coefficient
    Hd_B_            = squeeze(sum(Hd_B,2));
    Hd_A_            = squeeze(sum(Hd_A,2));
    Hd_B_            = Hd_B_';
    Hd_A_            = Hd_A_';
end


% calculate Power without Saturation
[power_sum,x_hat_F,power_F,F_PTO_F] = EqofMotion(Omega,Fex_,Mass,Hd_A_,...
    Hd_B_,Hs_K,PTO_damping,PTO_stiffness,pdata);


% Apply saturation conditions if necessary
if pdata.Opt.SatFlag

    if strcmpi(pdata.Opt.SatMethod,'PowerBound')

        % power in [W]
        power_sum = PowerBound(power_sum,pdata);

    end
end





% Review these before removing to avoid errors in all cases
%-------------------------------------------------------------------------%

% if  strcmpi(pdata.General.CaseStudy,'sim') && pdata.Opt.SatFlag
% 
%     pdata.Results.Force.F_PTO_sat2               = F_PTO_F_sat2;             % PTO force with Saturation
%     pdata.Results.Force.F_HS_sat2                = F_HS_F_sat2;              % Hydrostatic force with Saturation
%     pdata.Results.Force.F_Rad_sat2               = F_Rad_F_sat2;             % Radiation force with Saturation
%     pdata.Results.Power.P_sat                    = power_F_sat2;             % Power matrix with Saturation
%     pdata.Results.Power.power_sum_sat            = power_sum_sat2;           % Power sum for each device with Saturation
%     pdata.Results.Motion.x_sat                   = x_sat;                    % Heave displacement with Saturation
%     pdata.Results.Motion.u_sat                   = (1j.*Omega3).*x_sat;      % Heave velocity with Saturation
%     pdata.Results.Motion.a_sat                   = -(Omega3.^2).*x_sat;      % Heave acceleration with Saturation
%     pdata.Results.Motion.PTO_damping_sat_lim     = PTO_damping_sat_Lim;
%     pdata.Results.Motion.PTO_stiffness_sat_lim_L = PTO_stiffness_sat_Lim_L;
%     pdata.Results.Motion.PTO_stiffness_sat_lim_U = PTO_stiffness_sat_Lim_U;
%     power_sum                                    =  power_sum_sat2;
% end

% if strcmpi(pdata.General.CaseStudy,'sim')
%     % Hydrostatic and Hydrodynamic coefficients
%     pdata.Results.Coeffs.Hs_K    = Hs_K;       % Hydrostatic coeficient
%     pdata.Results.Coeffs.Hd_A    = Hd_A;       % Added mass
%     pdata.Results.Coeffs.Fe      = Fe;         % Excitation force
%     pdata.Results.Coeffs.Hd_B    = Hd_B;       % Damping coefficient
% 
%     % Forces
%     pdata.Results.Force.F_PTO      = F_PTO_F;         % PTO force
%     pdata.Results.Force.F_HS       = F_HS_F;          % Hydrostatic force
%     pdata.Results.Force.Fex_       = Fex_;            % Excitation force multiplied by amplitude
%     pdata.Results.Force.F_Rad      = F_Rad_F;         % Radiation force
% 
% 
%     % Power
%     pdata.Results.Power.P             = power_F;         % Power matrix
%     pdata.Results.Power.power_sum     = power_sum;       % Power sum for each device
% 
% 
%     % captor motion
%     pdata.Results.Motion.x       = x_hat_F;      % Heave displacement
%     pdata.Results.Motion.u       = u_hat_F;      % Heave velocity
%     pdata.Results.Motion.a       = a_hat_F;      % Heave acceleration
% 
% 
%     % Other terms
%     pdata.Results.other.Z_i      = Z_i;      % Intrinsic impedance
%     pdata.Results.other.Z_pto    = Z_pto_F;    % PTO impedance
%     pdata.Results.Mass = Mass;
% 
% end



end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [power_sum,x_hat_F,power_F, F_PTO_F,varargout] = EqofMotion(Omega,Fex_,...
    Mass,Hd_A_,Hd_B_,Hs_K,PTO_damping,PTO_stiffness,pdata)


% Ensure Omega is a row vector for power calculations
Omega = reshape(Omega,1,[]);

% Ensure Omega is consistent in size with the number of WECs
Omega = repmat(Omega,pdata.General.n_wec,1);

% Define intrinsic impedence
Z_i = Hd_B_ - (1j./Omega).*(-(Omega.^2).*(Mass + Hd_A_) + Hs_K);

% Check if Z_i is consistent in size
TF = isequal(size(Z_i), size(Fex_)) || (isvector(Z_i) &&....
    isvector(Fex_) && numel(Z_i) == numel(Fex_));

% Ensure Z_i3 is consistent in size 
if ~TF
    Z_i3 = repmat(Z_i,[1,1,length(pdata.WaveData.nodes(:,1))]);
else
    Z_i3 = Z_i;
end

% Define PTO impedence
Z_pto_F = PTO_damping - (1j./Omega).*PTO_stiffness;

% Check if Z_pto is consistent in size
TF = isequal(size(Z_pto_F), size(Fex_)) || (isvector(Z_pto_F) &&...
    isvector(Fex_) && numel(Z_pto_F) == numel(Fex_));

% Ensure Z_pto3 is consistent in size 
if ~TF
    Z_pto_F3 = repmat(Z_pto_F,[1,1,length(pdata.WaveData.nodes(:,1))]);
else
    Z_pto_F3 = Z_pto_F;
end

% Rearrange omega to avoid creating a for loop
Omega3 = repmat(Omega,1,1,length(pdata.WaveData.nodes(:,1)));

% Device motion in [m]
x_hat_F = (1./(1j.*Omega3)).*(Fex_./(Z_i3 + Z_pto_F3)); 
% x_hat_F_test = Fex_./(-(Omega3.^2)*(Mass+Hd_A_)-1j*Omega3*(PTO_damping+Hd_B_)+Hs_K);

% PTO force
[F_PTO_F,PTO_damping3] = CalcF_PTO(PTO_damping,PTO_stiffness,...
    Fex_, Omega3, x_hat_F, pdata);

% Power [w]
power_F = (1/2).*(PTO_damping3).*(Omega3.^2).*abs(x_hat_F).^2;
% power_F_text = (1/2).*(PTO_damping3).*(Omega3.^2).*abs(x_hat_F_test).^2;

if strcmpi(pdata.WaveData.Type,'irregular')
    
    % Add elements in the second dimension [W]
    DW = pdata.WaveData.DW;
    SP = pdata.WaveData.SP;
    power_sum = squeeze(sum(DW.*SP.*power_F,2));

    if pdata.General.n_wec == 1
        power_sum = power_sum';
    end

elseif strcmpi(pdata.WaveData.Type,'regular')
    % Add elements in the second dimension [W]
    power_sum  = power_F;                                                        
end


if nargout > 4

    % Device velocity in [m/s]: [u = (1j.*Omega).*x_hat]
    u_hat_F = (1j.*Omega3).*x_hat_F;

    % Device acceleration in [m/s^2]: [a = -(Omega.^2).*x_hat]
    a_hat_F          = -(Omega3.^2).*x_hat_F;

    % Hydrostatic force: [F_HS = -Hs_K.*x_hat]
    F_HS_F = -Hs_K.*x_hat_F;

    % Radiation force:
    Hd_B_3 = repmat(Hd_B_,[1,1,length(pdata.WaveData.nodes(:,1))]);
    Hd_A_3 = repmat(Hd_A_,[1,1,length(pdata.WaveData.nodes(:,1))]);
    F_Rad_F =  -1j.*Omega3.*Hd_B_3.*x_hat_F+ (Omega.^2).*Hd_A_3.*x_hat_F;

    varargout = {F_PTO_F,u_hat_F,a_hat_F,F_HS_F,F_Rad_F,Z_i,Z_pto_F,Omega3};

end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function [F_PTO_F,PTO_damping3,PTO_stiffness3] = CalcF_PTO(PTO_damping,...
    PTO_stiffness,Fex_,Omega3,x_hat_F,pdata)

% Check for size for PTO damping
TF = isequal(size(PTO_damping), size(Fex_)) || (isvector(PTO_damping)...
    && isvector(Fex_) && numel(PTO_damping) == numel(Fex_));

% Adjust size  for PTO damping
if ~TF
    PTO_damping3 = repmat(PTO_damping,[1,pdata.WaveData.nr,...
        length(pdata.WaveData.nodes(:,1))]);
else
    PTO_damping3 = PTO_damping;
end

% Check for size for PTO stiffness
TF = isequal(size(PTO_stiffness), size(Fex_))||(isvector(PTO_stiffness)...
    && isvector(Fex_) && numel(PTO_stiffness) == numel(Fex_));

% Adjust size  for PTO stiffness
if ~TF
    PTO_stiffness3   = repmat(PTO_stiffness,[1,pdata.WaveData.nr,...
        length(pdata.WaveData.nodes(:,1))]);
else
    PTO_stiffness3 = PTO_stiffness;
end

% PTO force calculations
F_PTO_F = -1j.*Omega3.*PTO_damping3.*x_hat_F - PTO_stiffness3.*x_hat_F;
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [power_sum, SatLim] = Saturations(Omega,Fex_,F_PTO_F,Hd_A_,...
    Hd_B_,Hs_K,x_hat_F,pdata)

PTO_damping = pdata.WEC.PTO_damping';
PTO_stiffness = pdata.WEC.PTO_stiffness';

% Plot: visualize device movement with no sturation
if strcmpi(pdata.General.WaveType,'irregular') && pdata.General.plotflag
    createSatPlots(Omega,F_PTO_F, x_hat_F, pdata, '0');
end


%------------------------ Stroke Limit -------------------------------- %
% WEC maximum displacement
x_max = pdata.WEC.Draft;

% Ratio of WEC displacement to its maximum
X_ratio = x_max./abs(x_hat_F);

% Saturation factor
rx = min(X_ratio,1);

% Ratio between the fundamental amplitude of the saturated
% waveform and the maximum displacement
alphax = (2/pi)*((1./rx).*asin(rx) + sqrt(1-(rx.^2)));

% Saturation factor
f_satx = alphax.*rx;

f_satx_min = min(f_satx(:));

if f_satx_min == 1 % Not saturated
    power_sum_sat = [];
    F_PTO_F_sat1  = [];
    x_sat         = x_hat_F;
else
    % Modify WEC displacement
    x_sat = f_satx.*x_hat_F;

    % Select one saturated PTO for each device using the worst case (minimum)?
    [power_sum_sat,F_PTO_F_sat1] =...
        CalcSatFX(Omega,x_sat,Fex_,Hd_A_, Hd_B_,Hs_K, PTO_damping,...
        PTO_stiffness,pdata);

    % Plot: visualize saturation considerations - device motion
    if strcmpi(pdata.General.WaveType,'irregular') && pdata.General.plotflag
        createSatPlots(Omega,F_PTO_F_sat1, x_sat, pdata, '2');
    end
end

%------------------------ PTO Force Limit -------------------------------- %

if isempty(F_PTO_F_sat1)
    F_PTO_F_sat1 = F_PTO_F;
end

% Ratio of PTO force to its maximum
F_ratio = abs(pdata.General.FMax./F_PTO_F_sat1);

% Saturation factor
rf = min(F_ratio,1);

% Ratio between the fundamental amplitude of the saturated
% waveform and the maximum force
alphaf = (2/pi)*((1./rf).*asin(rf) + sqrt(1-(rf.^2)));

% Saturation factor
f_satf = alphaf.*rf;

% find minimum saturation factor
f_satf_min = min(f_satf(:));

if f_satf_min == 1 % Not saturated

    for i =1:pdata.General.n_wec
        SatLim.PTO_damping_sat_Lim(i,1)     = pdata.Opt.B_PTO_max;
        SatLim.PTO_stiffness_sat_Lim_L(i,1) = pdata.Opt.k_PTO_min;
        SatLim.PTO_stiffness_sat_Lim_U(i,1) = pdata.Opt.k_PTO_max;
    end

    power_sum = power_sum_sat;

else

    % Modify Power
    PTO_damping_sat_     = f_satf.*PTO_damping;
    PTO_stiffness_sat_   = f_satf.*PTO_stiffness;

    %  Preallocate for faster computation
    PTO_damping_sat   = NaN(pdata.General.n_wec,1);
    PTO_stiffness_sat = NaN(pdata.General.n_wec,1);

    % Select one saturated PTO for each device using the worst case (minimum)?
    for i =1:pdata.General.n_wec
        PTO_damping_sat(i,1)     = min(min(PTO_damping_sat_(i,:,:)));
        PTO_stiffness_sat(i,1)   = min(min(abs(PTO_stiffness_sat_(i,:,:))))*sign(PTO_stiffness(i,1));
    end

    %  Preallocate for faster computation
    SatLim.PTO_damping_sat_Lim     = NaN(pdata.General.n_wec,1);
    SatLim.PTO_stiffness_sat_Lim_L = NaN(pdata.General.n_wec,1);
    SatLim.PTO_stiffness_sat_Lim_U = NaN(pdata.General.n_wec,1);

    % set upper limit to the bounds if not saturated
    for i =1:pdata.General.n_wec

        SatLim.PTO_damping_sat_Lim(i,1) = PTO_damping_sat(i,1);
        SatLim.PTO_stiffness_sat_Lim_L(i,1) = -abs(PTO_stiffness_sat(i,1));
        SatLim.PTO_stiffness_sat_Lim_U(i,1) = abs(PTO_stiffness_sat(i,1));

    end

    [power_sum_sat2,F_PTO_F_sat2] =...
        CalcSatFX(Omega,x_sat,Fex_,Hd_A_, Hd_B_,Hs_K, PTO_damping_sat,...
        PTO_stiffness_sat,pdata);

    if strcmpi(pdata.General.WaveType,'irregular') && pdata.General.plotflag
        createSatPlots(Omega,F_PTO_F_sat2, x_sat, pdata, '1');
    end

    power_sum = power_sum_sat2;
end






end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [power_sum,varargout] = CalcSatFX(Omega,x_sat,Fex_,Hd_A_, Hd_B_,...
    Hs_K, PTO_damping,PTO_stiffness,pdata)


% Ensure Omega is a row vector for power calculations
Omega = reshape(Omega,1,[]);

% 2-wec array Fex_ = [2 200]
Omega = repmat(Omega,pdata.General.n_wec,1);

% Rearrange omega to avoid creating a for loop
Omega3 = repmat(Omega,1,1,length(pdata.WaveData.nodes(:,1)));

[F_PTO_F,PTO_damping3] = CalcF_PTO(PTO_damping,...
    PTO_stiffness,Fex_,Omega3,x_sat,pdata);

% Hydrostatic force: [F_HS = -Hs_K.*x_hat]
F_HS_F = -Hs_K.*x_sat;   

% Radiation force 
Hd_B_3 = repmat(Hd_B_,[1,1,length(pdata.WaveData.nodes(:,1))]);
Hd_A_3 = repmat(Hd_A_,[1,1,length(pdata.WaveData.nodes(:,1))]);
F_Rad_F =  -1j.*Omega3.*Hd_B_3.*x_sat+ (Omega.^2).*Hd_A_3.*x_sat;

% Power in [W]
power_F = (1/2).*(PTO_damping3).*(Omega3.^2).*abs(x_sat).^2;

if strcmpi(pdata.WaveData.Type,'irregular')
    
    % Add elements in the second dimension [W]
    DW             = pdata.WaveData.DW;
    SP             = pdata.WaveData.SP;
    power_sum      = squeeze(sum(DW.*SP.*power_F,2));

    if pdata.General.n_wec == 1
        power_sum = power_sum';
    end

else
    % Add elements in the second dimension [W]
    power_sum  = sum(power_F);                                                        
end


if nargout == 2

    varargout = {F_PTO_F};

elseif nargout == 3

    varargout = {F_PTO_F,F_HS_F};

elseif nargout == 4

    varargout = {F_PTO_F,F_HS_F,F_Rad_F};

elseif nargout == 5

    varargout = {F_PTO_F,F_HS_F,F_Rad_F,power_F};

end


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function [power_sum] = PowerBound(power_sum,pdata)
% This function is created to apply an upper boud to the power that can 
% be generated by WECs. 
% The power upper bound is calculated in the CallWaveSpectrum function
% after the power elements are integrated over frequency domain. Therefore,
% this bound should beimposed not on the frequency-level elemnts of the
% power matrix, but the power elemnts at the sea state

if strcmpi(pdata.WaveData.Type,'regular')

    power_sum_ = power_sum;
    power_sum_(power_sum_>pdata.General.PowerLimit) = pdata.General.PowerLimit;

    % Add elements in the second dimension [W]
    power_sum  = power_sum_;

elseif strcmpi(pdata.WaveData.Type,'irregular')
    
    DW = pdata.WaveData.DW;
    SP = pdata.WaveData.SP;
    
    power_sum_ = power_sum;
    power_sum_(power_sum_>pdata.General.PowerLimit) = pdata.General.PowerLimit;


    if pdata.General.plotflag

        WEC_idx     = 1;

        % Modified (saturated) power
        pp_ = squeeze(power_sum_(WEC_idx,:));
        pp_m = reshape(pp_,[length(pdata.WaveData.nodes_Hs),length(pdata.WaveData.nodes_Te)]);

        % Power without saturation
        pp_0 = squeeze(power_sum(WEC_idx,:));
        pp_m0 = reshape(pp_0,[length(pdata.WaveData.nodes_Hs),length(pdata.WaveData.nodes_Te)]);

        CreateSurfPlot(pp_m,pp_m0, pdata)

    end

    % Power in [W]
    power_sum      = power_sum_;


    if pdata.General.n_wec == 1
        % Ensure row vector
        p = power_sum(:);
        power_sum = p';
    end

                                                     
end

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function createSatPlots(Omega,F_PTO_F, x_hat_F,pdata, Cflag)

commonFigureProperties;
switch Cflag
    case '0'
        hf = figure;
        CC = C.red(8,:);
    case '1'
        CC = C.blue(8,:);
        hf = gca;
    case '2'
        CC = C.green(8,:);
        hf = gca;
end



hf.Color = 'w';
subplot(1,2,1)
hold on;
for i=1:400
    plot(Omega,abs(F_PTO_F(1,:,i)),'.','Color',CC)
end
plot(Omega, pdata.General.FMax .*ones(size(Omega)),'--','Color', C.grey(10,:))
xlabel('$\omega$')
ylabel('$F_{pto}$')

subplot(1,2,2)
hold on;
for i=1:400
    plot(Omega,abs(x_hat_F(1,:,i)),'.','Color',CC)
end
plot(Omega, pdata.WEC.Draft.*ones(size(Omega)),'--','Color', C.grey(10,:))
xlabel('$\omega$')
ylabel('x')

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function CreateSurfPlot(pp_m,pp_m0,pdata)

[xq, yq]    = ndgrid(pdata.WaveData.nodes_Hs,pdata.WaveData.nodes_Te);

commonFigureProperties;
hf = figure;
hf.Color = 'w';

subplot(1,2,1)
surf(xq, yq,pp_m,'EdgeColor','none');
hold on;
l1 = length(pdata.WaveData.nodes_Hs);
l2 = length(pdata.WaveData.nodes_Te);
CO(:,:,1) = ones(l1,l2).*linspace(1,1,l2); % red
CO(:,:,2) = zeros(l1,l2).*linspace(0.5,0.6,l2); % green
CO(:,:,3) = zeros(l1,l2).*linspace(0,1,l2); % blue
surf(xq,yq,pdata.General.PowerLimit.*ones(size(pp_m)), CO,'EdgeColor','none','FaceAlpha',0.5)
xlabel('$H_s$', 'FontSize', fontlabel)
ylabel('$T_e$', 'FontSize', fontlabel)
zlabel('limited Power')
title('Saturated')
view(0,90)


subplot(1,2,2)
surf(xq, yq,pp_m0,'FaceAlpha',0.5,'EdgeColor','none')
hold on;
l1 = length(pdata.WaveData.nodes_Hs);
l2 = length(pdata.WaveData.nodes_Te);
CO(:,:,1) = ones(l1,l2).*linspace(1,1,l2); % red
CO(:,:,2) = zeros(l1,l2).*linspace(0.5,0.6,l2); % green
CO(:,:,3) = zeros(l1,l2).*linspace(0,1,l2); % blue
surf(xq,yq,pdata.General.PowerLimit.*ones(size(pp_m0)), CO,'EdgeColor','none','FaceAlpha',0.5)
xlabel('$H_s$', 'FontSize', fontlabel)
ylabel('$T_e$', 'FontSize', fontlabel)
zlabel('limited Power')
title('non-saturated')
view(0,90)

savename  = strcat(pdata.General.plot_dir, filesep, 'SaturationPlot','.pdf');
exportgraphics(hf, savename)
close(hf)

end






%-------------------------------------------------------------------------%
% Z_i_ = Hd_B - (1j./Omega).*(-(Omega.^2).*(pdata.Nemoh.Mass + Hd_A) + Hs_K);
% Z_pto_ = pdata.PTO.damping - (1j./Omega).*pdata.PTO.stiffness;
% x_ = (Fex_)./(pdata.PTO.stiffness + Hs_K - A_eff.*(Omega.^2) + 1j*(B_eff).*Omega);
% u_ = (1j.*Omega).*x_;
% a_ = -(Omega.^2).*x_;
% F_PTO_ = -1j*Omega*pdata.PTO.damping.*x_ - pdata.PTO.stiffness.*x_;
% F_HS_ = -Hs_K.*x_;
% F_Rad_ = -1j*Omega.*Hd_B.*x_ + (Omega.^2).*Hd_A.*x_;
% p_ = (1/2).*(pdata.PTO.damping).*(Omega.^2).*abs(x_).^2; 
%-------------------------------------------------------------------------%






















