function pdata = Process_Sim(pdata)

% This function takes as input certain conditions defined in Sim_Func to
% run various simulations. The current list in the switch statement below 
% corresponds to the list of simulations defined Sim_Func. 
% 
% 
% 
% Primary Contributer: Saeed Azad, PhD

switch pdata.SimCaseStudy

    case 'Case_1_Regular'

        pdata = Simulate_Regular(pdata);

    case 'Case_2_Regular'

         pdata = Simulate_Regular(pdata);

    case 'Case_Farm_landscape'

        pdata = Simulate_Farm_landscape(pdata);

    case 'Case_layout_sample'

        pdata = Simulate_Case_layout_sample(pdata);


    case 'Case_Load_Optim'

        % Make sure the correct power limit is used
        if pdata.Opt.SatFlag
            if pdata.sol.pdata.Opt.SatFlag == 1
                pdata.General.PowerLimit = pdata.sol.pdata.General.PowerLimit;
            end
            pdata.General.n_wec = 1;
            pdata = CallWaveSpectrum(pdata);
            pdata.WEC.x = pdata.WEC.x(1,1);
            pdata.WEC.y = pdata.WEC.y(1,1);
            pdata.WEC.PTO_damping = pdata.WEC.PTO_damping(1,1);
            pdata.WEC.PTO_stiffness = pdata.WEC.PTO_stiffness(1,1);
            pdata.WEC.Distance = 0;
            pdata.WEC.Angle = 0;
        end
        
        % Dummy variable
        pdata0 = pdata;


        % Run simulation using MS
        pdata0.HydroFlag = 'MS';
        MS_sim           = run_simulation(pdata0);


        % Average power per WEC per year ub kW
        MS_sim.Results.Avg_Power_wec_year_kW_MS = ((sum(MS_sim.Results.Power.P_total_vec,2))/1000)./30;

        % Calculate q-factor
        [Q_factor, q_factor]               = calculate_q_factor(MS_sim);
        MS_sim.Results.Q_factor_MS         = Q_factor;
        MS_sim.Results.q_factor_MS         = q_factor;

        pdata.Results.MS                   = MS_sim.Results;

        savename = strcat(pdata.General.Solution_dir,filesep,'Solution_',pdata.HydroFlag,'_',pdata0.WaveData.Region);
        save(savename,'pdata')

    case 'Case_Radius_v_kPTO'

        % run simulation using surrogate models in MBE
        pdata = Case_Radius_v_kPTO(pdata);

    otherwise


        % Run simulation using MS
        if strcmpi(pdata.HydroFlag,'MS')

            Sol = run_simulation(pdata);
            
            % Calculate q-factor
            [Q_factor, q_factor]            = calculate_q_factor(Sol);
            Sol.Results.Q_factor_MS         = Q_factor;
            Sol.Results.q_factor_MS         = q_factor;

        elseif strcmpi(pdata.HydroFlag,'SM_MBE')

            % Use SM_MBE
            SM_MBE = run_simulation(pdata);

            % Calculate q-factor
            [Q_factor, q_factor]               = calculate_q_factor(SM_MBE);
            SM_MBE.Results.Q_factor_MS         = Q_factor;
            SM_MBE.Results.q_factor_MS         = q_factor;

            % Use MS
            pdata0 = pdata;
            pdata.HydroFlag = 'MS';
            MS = run_simulation(pdata);
            
            % Calculate q-factor
            [Q_factor, q_factor]           = calculate_q_factor(MS);
            MS.Results.Q_factor_MS         = Q_factor;
            MS.Results.q_factor_MS         = q_factor;

            Sol.SM_MBE = SM_MBE;
            Sol.MS = MS;

            pdata = pdata0;
        end


        savename = strcat(pdata.General.Solution_dir,filesep,filesep,'Solution_',pdata.HydroFlag,'_',pdata.WaveData.Region);
        save(savename,'Sol')
        

end

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = Simulate_Regular(pdata)

% run simulation
pdata        = run_simulation(pdata);

% Power [kW]
Power         = sum(sum(pdata.Results.Power.P_total_vec))/1000;    
pdata.Results.Power.Power_total = Power;

savename = strcat(pdata.General.Solution_dir,filesep,'Sim_sol');
save(savename,'pdata');

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function pdata = Simulate_Farm_landscape(pdata)

% Define dummy variable
X   = pdata.WEC.x;
Y   = pdata.WEC.y;

% Dummy variable - needed for parallel implementation
P_DATA = pdata;

% number of simulations
n_sim  = pdata.n_sim;
dis    = pdata.WEC.Distance;
ang    = pdata.WEC.Angle;

% Pre-assign
Power  = NaN(1,n_sim);
Q_fac  = cell(1,n_sim);
q_fac  = cell(1,n_sim);


% Perform simulation
parfor i = 1:n_sim

    % Dummy variable - broadcast
    P_DATA               = pdata;
    P_DATA.WEC.x         = X(:,i);
    P_DATA.WEC.y         = Y(:,i);
    P_DATA.WEC.Distance  = dis(:,i);
    P_DATA.WEC.Angle     = ang(:,i);

    % run simulation
    Sim_sol         = run_simulation(P_DATA);

    % Energy [kW]
    Power(1,i)     = sum(sum(Sim_sol.Results.Power.P_total_vec))/1000;


    % Calculate q-factor
    [Q_factor, q_factor]     = calculate_q_factor(Sim_sol);

    Q_fac{1,i} = Q_factor;
    q_fac{:,i} = q_factor;

end

Q.Q_fac = Q_fac;
Q.q_fac = q_fac;

x_vec_p                 = pdata.WEC.x_vec_p;
y_vec_p                 = pdata.WEC.y_vec_p;
x_vec_n                 = -x_vec_p;
y_vec_n                 = -y_vec_p;
x_vec                   = [fliplr(x_vec_n),x_vec_p];
y_vec                   = [fliplr(y_vec_n),y_vec_p];
[x_mat,y_mat]           = meshgrid(x_vec,y_vec);

savenameE = strcat(pdata.General.Solution_dir,filesep,'DS_Power_',pdata.HydroFlag,string(max(x_vec_p)));
save(savenameE,'Power');

savenameQ = strcat(pdata.General.Solution_dir,filesep,'DS_Q_',pdata.HydroFlag,string(max(x_vec_p)));
save(savenameQ,'Q');


% Create Plots - q-factor
q_ = Q.Q_fac;
q  = cell2mat(q_);
q_mat = reshape(q,2*length(pdata.WEC.y_vec_p),2*length(pdata.WEC.x_vec_p));
contourf(x_mat ,y_mat,q_mat,20,'edgecolor','none')
colormap(parula)
c = colorbar;
clim([0.9, 1.05])
xlabel('x')
ylabel('y')
c.Label.String = 'q_{factor}';
savefig(strcat(pdata.General.Solution_dir,filesep,'Q_R5',string(max(x_vec_p))))
exportgraphics(gca,strcat(pdata.General.Solution_dir,filesep,'Q_R5',string(max(x_vec_p)),'.pdf'));

% Create Plots - energy
figure(2)
e_mat = reshape(Power,length(y_vec),length(x_vec));
contourf(x_mat,y_mat,e_mat,20,'edgecolor','none')
colormap(parula)
c = colorbar;
xlabel('x')
ylabel('y')
c.Label.String = 'Power [kW]';
clim([130, 200])
savefig(strcat(pdata.General.Solution_dir,filesep,'Energy_',string(max(x_vec_p))))
exportgraphics(gca,strcat(pdata.General.Solution_dir,filesep,'Energy_',string(max(x_vec_p)),'.pdf'))




end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = Simulate_Case_layout_sample(pdata)

% Layout variables
X   = pdata.WEC.x;
Y   = pdata.WEC.y;
Dis = pdata.WEC.Distance;
Ang = pdata.WEC.Angle;

% Plant variables
R  = pdata.WEC.Radius;
RD = pdata.WEC.RD;
Draft = pdata.WEC.Draft;

% Control variables
B_pto = pdata.WEC.PTO_damping;
K_pto = pdata.WEC.PTO_stiffness;


% Dummy variable - needed for parallel implementation
P_DATA = pdata;

% number of simulations
N_l  = pdata.N_l;


% Pre-assign
Power = NaN(1,N_l);
Q_fac  = cell(1,N_l);
q_fac  = cell(1,N_l);
    

t1 = tic;

% Perform simulation
for i = 1:N_l

    % Broadcast variable
    P_DATA               = pdata;

    % Layout variables 
    P_DATA.WEC.x         = X(:,i);
    P_DATA.WEC.y         = Y(:,i);
    P_DATA.WEC.Distance  = Dis(:,i);
    P_DATA.WEC.Angle     = Ang(:,i);

    % Plant
    P_DATA.WEC.Radius    = R(1,i);
    P_DATA.WEC.RD        = RD(1,i);
    P_DATA.WEC.Draft     = Draft(1,i);

    % Control
    P_DATA.WEC.PTO_damping = B_pto(i,:);
    P_DATA.WEC.PTO_stiffness = K_pto(i,:);


    % run simulation
    Sim_sol         = run_simulation(P_DATA);

    % Power [W]
    Power(1,i)     = sum(sum(Sim_sol.Results.Power.P_total_vec));


    % Calculate q-factor
    [Q_factor, q_factor]     = calculate_q_factor(Sim_sol);

    Q_fac{1,i} = Q_factor;
    q_fac{:,i} = q_factor;

end

Q.Q_fac = Q_fac;
Q.q_fac = q_fac;

Volume = pi*(R.^2)*2.*Draft;
Obj    = Power/Volume; % W/m^3

t2 = toc(t1);


switch pdata.General.run_case

    case 'batch_simulations'

        savenameE = strcat(pdata.General.Solution_dir,filesep,'Power_',pdata.HydroFlag,'_',string(pdata.N_l));
        save(savenameE,'Power');

        savenameQ = strcat(pdata.General.Solution_dir,filesep,'Q_',pdata.HydroFlag,'_',string(pdata.N_l));
        save(savenameQ,'Q');

        savenameObj = strcat(pdata.General.Solution_dir,filesep,'Obj_',pdata.HydroFlag,'_',string(pdata.N_l));
        save(savenameObj,'Obj');


    case 'run_optimized'
        savenameE = strcat(pdata.General.Solution_dir,filesep,'power_',pdata.HydroFlag,'_',string(pdata.N_l),'_',pdata.General.run_case);
        save(savenameE,'Power');

        savenameQ = strcat(pdata.General.Solution_dir,filesep,'Q_',pdata.HydroFlag,'_',string(pdata.N_l),'_',pdata.General.run_case);
        save(savenameQ,'Q');

end
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function pdata = Case_Radius_v_kPTO(pdata)

% Make sure the correct power limit is used
if pdata.Opt.SatFlag
    pdata.General.PowerLimit = pdata.General.PowerLimit;
end

% Layout variables
X   = pdata.WEC.x;
Y   = pdata.WEC.y;


% Plant variables
R     = pdata.WEC.Radius;
RD    = pdata.WEC.RD;
Draft = pdata.WEC.Draft;

% Control variables
B_pto = pdata.WEC.PTO_damping;
K_pto = pdata.WEC.PTO_stiffness;


% Dummy variable - needed for parallel implementation
P_DATA = pdata;

% number of simulations
N_l  = pdata.n_sim;


% Pre-assign
Obj    = NaN(N_l,1);
Power  = NaN(N_l,1);

Volume =  pdata.General.n_wec*(pi*R.^2)*2.*Draft;

t1 = tic;

% Perform simulation
% parfor (i = 1:N_l,40)
for i = 1:N_l

    % Broadcast variable
    P_DATA               = pdata;

    % Layout variables 
    P_DATA.WEC.x         = X;
    P_DATA.WEC.y         = Y;

    % Plant
    P_DATA.WEC.Radius    = R(i,1);
    P_DATA.WEC.RD        = RD;
    P_DATA.WEC.Draft     = Draft(i,1);

    % Control
    P_DATA.WEC.PTO_damping = B_pto(i,:);
    P_DATA.WEC.PTO_stiffness = K_pto(i,:);


    % run simulation
    Sim_sol         = run_simulation(P_DATA);

    % Power [W]
    Power(i,1)     = sum(sum(Sim_sol.Results.Power.P_total_vec))/1000;

    % Objective funciton
    Obj(i,1)       = Power(i,1)/Volume(i,1);


end

t2 = toc(t1);


R_matrix = reshape(R,pdata.General.n1,pdata.General.n2,pdata.General.n3);
k_matrix = reshape(K_pto(:,1),pdata.General.n1,pdata.General.n2,pdata.General.n3);
B_matrix = reshape(B_pto(:,1),pdata.General.n1,pdata.General.n2,pdata.General.n3);
Obj_matrix  = reshape(Obj,pdata.General.n1,pdata.General.n2,pdata.General.n3);
power_matrix = reshape(Power,pdata.General.n1,pdata.General.n2,pdata.General.n3);

Sim_results.R_matrix = R_matrix;
Sim_results.k_matrix = k_matrix;
Sim_results.B_matrix = B_matrix;
Sim_results.Obj_matrix = Obj_matrix;
Sim_results.power_matrix = power_matrix;


% Plot results: radius versus B_pto
x = squeeze(B_matrix(:));
y = squeeze(k_matrix(:));
z = squeeze(R_matrix(:));
o = squeeze(Obj_matrix(:));

commonFigureProperties;
hf = figure;
hf.Color = 'w';
scatter3(x,y,z,1,o,'filled');
ax = gca;
% ax.XDir = 'reverse';
view(-31,14)
xlabel('$B_{pto}$')
ylabel('$K_{pto}$')
zlabel('Radius')
colormap(turbo)
cb = colorbar;
cb.Label.String = 'Objective';
hold on

for  i= 1:pdata.General.n1
     mat_ = squeeze(Obj_matrix(i,:,:));
     maximum = max(max(mat_));
    [r,c]    = find(mat_==maximum);
    B_pto_opt(i,1) = B_matrix(i,r,c);
    K_pto_opt(i,1) = k_matrix(i,r,c);
    R_opt(i,1)     = R_matrix(i,1,1);
    Obj_opt(i,1)   = maximum;
end

n_ = 1;
for  i= 1:pdata.General.n1-1
    rl(i,:) = linspace(R_opt(i,1),R_opt(i+1),n_);
    Bl(i,:) = linspace(B_pto_opt(i,1),B_pto_opt(i+1),n_);
    Kl(i,:) = linspace(K_pto_opt(i,1),K_pto_opt(i+1),n_);
end

scatter3(Bl,Kl,rl,'*','MarkerEdgeColor',C.purple(5,:),'MarkerFaceColor',C.lime(5,:))


if pdata.Opt.SatFlag
    nn = string(pdata.General.PowerLimit); 
else
    nn = 'none'; 
end

savenameE = strcat(pdata.General.Solution_dir,filesep,'Sim_results_',nn);
save(savenameE,'Sim_results');

savename_ = strcat(pdata.General.Solution_dir,filesep,'RvK_Sim_results_',nn,'.pdf');
exportgraphics(gca,savename_)

savename_ = strcat(pdata.General.Solution_dir,filesep,'RvK_Sim_results_',nn,'.fig');
saveas(gca,savename_)


end
