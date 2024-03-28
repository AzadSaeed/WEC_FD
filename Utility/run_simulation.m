function pdata = run_simulation(pdata,varargin)
% This function run various simulations of WEC farms based on
% specifications defined in pdata. The second input argument is a flag that
% when is set to 'Q_factor' calculates the q-factor 

if nargin == 1
    Flag = [];
elseif nargin ==2
    Flag = varargin{1,1};
end


% Run simulation
SimulationTime = tic;
[Obj, pdata, P_total_vec, power_, p_mat, Pe, Pae] = OBJ_Calculations([],pdata);
pdata.Results.Sim_Time = toc(SimulationTime);


% Store simulation results

 % Power for all sea states
pdata.Results.Power.power_ = power_;          

% Power matrix without climate probability or limits
pdata.Results.Power.p_mat = p_mat;       

% Power matrix accounting for PCC efficiency
pdata.Results.Power.Pe = Pe;

% probabilistic power for all years
pdata.Results.Power.Pae = Pae;     

% Total yearly power with all efficiency factors in [W]
pdata.Results.Power.P_total_vec = P_total_vec;     


% Objective function [kW/m^3]
if ~strcmpi(pdata.General.CaseStudy,'Sim')
    pdata.Results.Power.Obj = Obj;
end


% Organize Outputs in a Table
if ~strcmpi(pdata.General.CaseStudy,'Sim')
    if ~strcmpi(Flag,'Q_factor') || isempty(Flag)
        T = TableReport(pdata);
        pdata.Results.T = T;
    end
end




% if pdata.General.plotflag
%     switch pdata.General.CaseStudy
% 
%         case 'PCL_opt'
%             % Create_Sim_plot(pdata)
%             close all
%             Create_Ary_plot(pdata)
%             Create_Ary_Cons(pdata)
% 
%         case 'C_opt'
%             close all
%             Create_Ary_plot(pdata)
% 
%         case 'A_opt'
%             close all
%             Create_Ary_plot(pdata)
%             Create_Ary_Cons(pdata)
% 
%         case 'Sim'
% 
%             if strcmpi(pdata.General.WaveType,'Regular')
% 
%                 Regular_Time_Plots(pdata);
% 
%             end
% 
%     end
%     close all;
% end

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function Create_Ary_plot(pdata)


commonFigureProperties;
% hf = figure;
% hf.Color = 'w';
% hold on;

D = pdata.WEC.Distance;
theta = pdata.WEC.Angle;
polarplot(theta,D,'o','LineWidth',linewidth, 'MarkerSize', 10, 'MarkerEdgeColor', C.blue(9,:),'MarkerFaceColor', C.yellow(8,:))

savename = strcat(pdata.General.plot_dir, filesep,'OptimalArray_',string(pdata.Opt.MaxFuncEval),'.pdf');
exportgraphics(gca,savename);
close all;

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function Create_Ary_Cons(pdata)


commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on;

Sd = pdata.General.SafeDist;

rr = (2*pdata.WEC.Radius + Sd)/2;
XL = [pdata.WEC.x, pdata.WEC.y];

for jj = 1:length(XL)
    viscircles(XL(jj,:),rr,'LineStyle','--', 'color', 'r')
end
hold on
plot(XL(:,1),XL(:,2),'o','LineWidth',linewidth,'MarkerEdgeColor', C.blue(9,:),'MarkerFaceColor', C.yellow(8,:))
% area = 0.5*sqrt(20000*pdata.General.n_wec);
rectangle('Position',[pdata.Opt.lbl(1,1) pdata.Opt.lbl(2,1) pdata.Opt.ubl(1,1) 2*pdata.Opt.ubl(1,1)],'LineStyle','-.')
xlim([-300 300])
ylim([-300 300])
axis equal
savename = strcat(pdata.General.plot_dir, filesep,'OptimalArrayConstraint',string(pdata.General.n_wec),'_',string(pdata.Opt.MaxFuncEval),'.pdf');
exportgraphics(gca,savename);
close(hf)

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Regular_Time_Plots(pdata)

Hs                  = pdata.WaveData.Hs_regular;
Te                  = pdata.WaveData.Te_regular;
Omega               = (2*pi)/Te;
t                   = linspace(0,60,3600);

x_hat               = pdata.Results.Motion.x;
u_hat               = pdata.Results.Motion.u;
a_hat               = pdata.Results.Motion.a;

xt_                 = real(x_hat.*exp(1j.*Omega.*t));
ut_                 = real(u_hat.*exp(1j.*Omega.*t));
at_                 = real(a_hat.*exp(1j.*Omega.*t));

pdata.Results.Time_series.Omega = Omega;
pdata.Results.Time_series.t     = t;
pdata.Results.Time_series.xt    = xt_;
pdata.Results.Time_series.ut    = ut_;
pdata.Results.Time_series.at    = at_;


% Create time-domain plots
Createtimeplots(pdata, 'x')
Createtimeplots(pdata, 'v')
Createtimeplots(pdata, 'a')
close all

if pdata.Opt.SatFlag
    x_hat_sat           = pdata.Results.Motion.x_sat;
    u_hat_sat           = pdata.Results.Motion.u_sat;
    a_hat_sat           = pdata.Results.Motion.a_sat;

    xt_sat              = real(x_hat_sat.*exp(1j.*Omega.*t));
    ut_sat              = real(u_hat_sat.*exp(1j.*Omega.*t));
    at_sat              = real(a_hat_sat.*exp(1j.*Omega.*t));

    pdata.Results.Time_series.Omega = Omega;
    pdata.Results.Time_series.t     = t;
    pdata.Results.Time_series.xt    = xt_sat;
    pdata.Results.Time_series.ut    = ut_sat;
    pdata.Results.Time_series.at    = at_sat;

    Createtimeplots(pdata, 'x')
    Createtimeplots(pdata, 'v')
    Createtimeplots(pdata, 'a')
end

% Create mean power input versus k and versus B
CreatePvKBplots(pdata,'pk')
CreatePvKBplots(pdata,'pB')
close all

% Plot Regions of Resonance
CreateRRplots(pdata,'Prescribed')
CreateRRplots(pdata,'Optimal')



end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Createtimeplots(pdata, plotflag)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on;

idx = 1; % desired node
t   = pdata.Results.Time_series.t(1,:,idx);

switch plotflag
    
    case 'x'

        switch upper(pdata.WaveData.Type)
            case 'REGULAR'
                x = pdata.Results.Time_series.xt;
            case 'IRREGULAR'
                x = pdata.Results.Time_series.xt(1,:,idx);
        end
        plot(t,x,'LineWidth',linewidth, 'MarkerSize', 12, 'Color', C.red(9,:))
        titletext = 'Buoy displacement [$m$]';
        yl = '$x~[m]$';

    case 'v'

        switch upper(pdata.WaveData.Type)
            case 'REGULAR'
                x = pdata.Results.Time_series.ut;
            case 'IRREGULAR'
                x = pdata.Results.Time_series.ut(1,:,idx);
        end
        plot(t,x,'LineWidth',linewidth, 'MarkerSize', 12, 'Color', C.blue(9,:))
        titletext = 'Buoy velocity [$m/s$]';
        yl = '$v~[m/s]$';

    case 'a'

        switch upper(pdata.WaveData.Type)
            case 'REGULAR'
                x = pdata.Results.Time_series.at;
            case 'IRREGULAR'
                x = pdata.Results.Time_series.at(1,:,idx);
        end
        plot(t,x,'LineWidth',linewidth, 'MarkerSize', 12, 'Color', C.green(9,:))
        titletext = 'Buoy acceleration [$m/s^2$]';
        yl = '$a~[m/s^2]$';

end

ax = gca;
ax.FontSize = fonttick;
xlabel('$t~[s]$', 'FontSize', fontlabel)
ylabel(yl, 'FontSize', fontlabel)
title(titletext)
savename = strcat(pdata.General.s_path, filesep, strcat(plotflag,'.pdf'));
taskflag = 'axes';
commonFigureTasks;
taskflag = 'legend';
commonFigureTasks;
exportgraphics(gca,savename);


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function CreateFrplots(pdata, plotflag)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on;

idx = 1; % desired node

switch plotflag

    case 'x'

        w = pdata.Results.Time_series.Omega;
        Z = pdata.Results.Motion.x;
        plot(w,abs(Z(:,idx)),'.-','LineWidth',linewidth, 'MarkerSize', 12, 'Color', C.red(9,:))
        yl = '$\vert x \vert$';

    case 'v'

        w = pdata.Results.Time_series.Omega;
        Z_dot = pdata.Results.Motion.u;
        plot(w,abs(Z_dot(:,idx)),'.-','LineWidth',linewidth, 'MarkerSize', 12, 'Color', C.blue(9,:))
        yl = '$ \vert \dot{x} \vert$';

    case 'a'

        w = pdata.Results.Time_series.Omega; 
        Z_ddot = pdata.Results.Motion.a; 
        plot(w,abs(Z_ddot(:,idx)),'.-','LineWidth',linewidth, 'MarkerSize', 12, 'Color', C.green(9,:))
        yl = '$ \vert \ddot{x} \vert$';

end

ax = gca;
ax.FontSize = fonttick;
xlabel('$\omega$', 'FontSize', fontlabel)
ylabel(yl, 'FontSize', fontlabel)
savename = strcat(pdata.General.s_path, filesep, strcat('f',plotflag,'.pdf'));
taskflag = 'axes';
commonFigureTasks;
taskflag = 'legend';
commonFigureTasks;
exportgraphics(gca,savename);


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function  CreatePvKBplots(pdata, plotflag)
% This function plots power against a range of k and B values
% It is mainly developed for frequency-domain analysis with regular waves


Hs                  = pdata.WaveData.Hs_regular;
Te                  = pdata.WaveData.Te_regular;
Omega               = (2*pi)/Te;
Fex_                = pdata.Results.Force.Fex_;


Z_i    = pdata.Results.other.Z_i;
Z_pto  = @(k,B) B - (1j/Omega)*k;
x_hat  = @(k,B) (1/(1j*Omega)).*(Fex_./(Z_i + Z_pto(k,B)));
u_hat  = @(k,B) (1j.*Omega).*x_hat(k,B);
power  = @(k,B) (1/2).*(B).*(Omega.^2).*abs(x_hat(k,B)).^2;
commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on;


switch plotflag

    case 'pk'

        k                = linspace(-2*10^6,10^6,10000);                             % Range of k values
        power_val_p      = power(k,pdata.WEC.PTO_damping)/1000;                        % Power in [kN.(m/S).(rad^2)]
        power_val_opt    = power(k,pdata.Results.PTO.B_optimal)/1000;              % Power in [kN.(m/S).(rad^2)]
        Power_opt        = power(pdata.Results.PTO.k_optimal,pdata.Results.PTO.B_optimal)/1000;
        x                = k;
        x3               = pdata.Results.PTO.k_optimal;
        leg              = {'Power with prescribed $B_{PTO}$',...
                            'Power with optimal $B_{PTO}$',...
                            'Power with optimal $k_{PTO}$ and $B_{PTO}$'};
        titletext        = 'Power versus $k_{pto}$';
        xL               = '$k_{pto}~[N/m]$';

    case 'pB'

        B                = linspace(0,10^6,10000);                                 % Range of B values 
        power_val_p      = power(pdata.WEC.PTO_stiffness,B)/1000;                      % Power in [kN.(m/S).(rad^2)]
        power_val_opt    = power(pdata.Results.PTO.k_optimal,B)/1000;              % Power in [kN.(m/S).(rad^2)]
        Power_opt        = power(pdata.Results.PTO.k_optimal,pdata.Results.PTO.B_optimal)/1000;
        x                = B;
        x3               = pdata.Results.PTO.B_optimal;
        leg              = {'Power with prescribed $k_{PTO}$',...
                            'Power with optimal $k_{PTO}$',...
                            'Power with optimal $k_{PTO}$ and $B_{PTO}$'};
        titletext        = 'Power versus $B_{pto}$';
        xL               = '$B_{pto}~[N.s/m]$';

end



plot(x,power_val_p,'LineWidth',linewidth, 'MarkerSize', 12,...
    'Color', C.blue(7,:))
plot(x,power_val_opt,'LineWidth',linewidth, 'MarkerSize', 12,...
    'Color', C.red(7,:))
plot(x3,Power_opt,'*','LineWidth',linewidth,...
    'MarkerSize', 12, 'Color', C.lightgreen(9,:),'MarkerEdgeColor',...
    C.lightgreen(8,:), 'MarkerFaceColor',C.lightgreen(8,:) )


%title(titletext)
legend(leg,'Location','northoutside')        
xlabel(xL, 'FontSize', fontlabel)
ylabel('$Power~[kW]$', 'FontSize', fontlabel)

ax = gca;
ax.FontSize = fonttick;
savename = strcat(pdata.General.s_path, filesep, strcat(plotflag,'.pdf'));

taskflag = 'axes';
commonFigureTasks;
taskflag = 'legend';
commonFigureTasks;
exportgraphics(gca,savename);


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function CreateAEPplots(AEP_vec, pdata)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on;
plot(1:length(AEP_vec), AEP_vec/1000,'.-',...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.red(9,:),'MarkerEdgeColor',C.blue(8,:))
ax = gca;
ax.FontSize = fonttick;
xlabel('$Year$', 'FontSize', fontlabel)
ylabel('$Energy$', 'FontSize', fontlabel)
title('Annual Energy Production [MWh]')
savename = strcat(pdata.General.s_path, filesep, strcat('AEP','.pdf'));
taskflag = 'axes';
commonFigureTasks;
taskflag = 'legend';
commonFigureTasks;
exportgraphics(gca,savename);

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function CreateRRplots(pdata,plotflag)
% This function plots regions of resonace.
% For the interpretation of this plot refer to 
% https://www.acs.psu.edu/drussell/Demos/Resonance-Regions/Resonance.html

commonFigureProperties;
hf = figure;
hf.Color = 'w';


w         = linspace(0.0001,3,1000);
Fex_      = pdata.Results.Force.Fex_;
Hd_A      = pdata.Results.Coeffs.Hd_A;
Hd_B      = pdata.Results.Coeffs.Hd_B;
Hs_K      = pdata.Results.Coeffs.Hs_K;
Z_i       = @(W) Hd_B - (1j./W).*(-(W.^2).*(pdata.Results.Mass + Hd_A) + Hs_K);     % Intrinsic Impedance: [Z_i = Hd_B - (1j./Omega).*(-(Omega.^2).*(pdata.Nemoh.Mass + Hd_A) + Hs_K)]
Z_pto     = @(W,k,B) B - (1j./W)*k;                                               % PTO impedance: [Z_pto = pdata.PTO.damping - (1j./Omega).*pdata.PTO.stiffness]
x_hat     = @(W,k,B) (1./(1j*W)).*(Fex_./(Z_i(W) + Z_pto(W,k,B)));                % Captor motion: [x = (a_i.*Fex)./(K - A.*(Omega.^2) + 1j*B.*Omega)];
power     = @(W,k,B) (1/2).*(B).*(W.^2).*abs(x_hat(W,k,B)).^2; 

switch plotflag

    case'Prescribed'

        pdata.Results.PTO_stiffness = pdata.Results.PTO.k_prescribed;
        pdata.Results.Radius = pdata.WEC.Radius;
        pdata.Results.Draft  = pdata.WEC.Draft;
        Omega_n = Calc_Natural_Frequency(pdata);
        pdata.Results.Omega_n  = Omega_n;
        omega    = w./pdata.Results.Omega_n;
        x        = abs(x_hat(omega,pdata.WEC.PTO_stiffness,pdata.WEC.PTO_damping));
        p        = power(omega,pdata.WEC.PTO_stiffness,pdata.WEC.PTO_damping);
        savename = strcat(pdata.General.s_path, filesep, strcat('RofR_prs',plotflag,'.pdf'));
        B_test   = 10^3;
        K_test   = -10^5;
        x_Btest  = abs(x_hat(omega,pdata.WEC.PTO_stiffness,B_test));
        pB_test  = power(omega,pdata.WEC.PTO_stiffness,B_test);
        x_ktest  = abs(x_hat(omega,K_test,pdata.WEC.PTO_damping));
        pk_test  = power(omega,K_test,pdata.WEC.PTO_damping);
        leg_text1 = sprintf('$k_{pto} = %.1e$ and $B_{pto} = %.1e$',pdata.WEC.PTO_stiffness,pdata.WEC.PTO_damping);
        leg_text2 = sprintf('$k_{pto} = %.1e$ and $B_{pto} = %.1e$',pdata.WEC.PTO_stiffness,B_test);
        leg_text3 = sprintf('$k_{pto} = %.1e$ and $B_{pto} = %.1e$',K_test,pdata.WEC.PTO_damping);
        leg_text  = {leg_text1,leg_text2,leg_text3};
    case 'Optimal'
        pdata.Results.PTO_stiffness = pdata.Results.PTO.k_optimal;
        pdata.Results.Radius = pdata.WEC.Radius;
        pdata.Results.Draft  = pdata.WEC.Draft;
        Omega_n = Calc_Natural_Frequency(pdata);
        pdata.Results.Omega_n  = Omega_n;
        omega    = w./pdata.Results.Omega_n;
        x        = abs(x_hat(omega,pdata.Results.PTO.k_optimal,pdata.Results.PTO.B_optimal));
        p        = power(omega,pdata.Results.PTO.k_optimal,pdata.Results.PTO.B_optimal);
        savename = strcat(pdata.General.s_path, filesep, strcat('RofR_opt',plotflag,'.pdf'));
        B_test   = 10^6;
        K_test   = 10^6;
        x_Btest  = abs(x_hat(omega,pdata.Results.PTO.k_optimal,B_test));
        pB_test  = power(omega,pdata.Results.PTO.k_optimal,B_test);
        x_ktest  = abs(x_hat(omega,K_test,pdata.Results.PTO.B_optimal));
        pk_test  = power(omega,K_test,pdata.Results.PTO.B_optimal);
        leg_text1 = sprintf('$k_{pto} = %.1e$ and $B_{pto} = %.1e$',pdata.Results.PTO.k_optimal,pdata.Results.PTO.B_optimal);
        leg_text2 = sprintf('$k_{pto} = %.1e$ and $B_{pto} = %.1e$',pdata.Results.PTO.k_optimal,B_test);
        leg_text3 = sprintf('$k_{pto} = %.1e$ and $B_{pto} = %.1e$',K_test,pdata.Results.PTO.B_optimal);
        leg_text  = {leg_text1,leg_text2,leg_text3};
end

semilogy(omega,x,...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.blue(9,:),'MarkerEdgeColor',C.blue(8,:))
hold on
semilogy(omega,x_Btest,...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.blue(5,:),'MarkerEdgeColor',C.blue(5,:))

semilogy(omega,x_ktest,...
    'linewidth',linewidth,'markersize',12,...
    'Color',C.blue(3,:),'MarkerEdgeColor',C.blue(3,:))

ax = gca;
ax.FontSize = fonttick;
xlabel('$\frac{\omega}{\omega_{n}}$', 'FontSize', fontlabel)
ylabel('$\vert x \vert$', 'FontSize', fontlabel)
legend(leg_text,'Location','northoutside')
taskflag = 'legend';
commonFigureTasks;
exportgraphics(gca,savename);


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function T = TableReport(pdata)

% Fields = ['Radius', 'Slenderness', 'Draft', 'PTO Damping', 'PTO Stiffness', 'Power', 'PowerPerVolume'];

Radius = pdata.Results.Radius*ones(pdata.General.n_wec,1);
Slenderness = pdata.Results.RD*ones(pdata.General.n_wec,1);
Draft = pdata.Results.Draft*ones(pdata.General.n_wec,1);
PTO_damping = pdata.Results.PTO_damping';
PTO_stiffness = pdata.Results.PTO_stiffness';
Power_vec = sum(pdata.Results.Power.P_total_vec,2);
Volume = pi*(pdata.Results.Radius^2)*Draft*2;
PowerPerVolume = Power_vec./Volume;

T = table(Radius, Slenderness, Draft, PTO_damping, PTO_stiffness, Power_vec, PowerPerVolume);


end

