function [Q_factor, q_factor] = calculate_q_factor(pdata)

% Create Dummy varibale
pd                  = pdata;

% Set the numbers of WECs to one
pd.General.n_wec    = 1;

% calculate the wave spectrum again 
pd                  = CallWaveSpectrum(pd);

pd.WEC.Radius       = pdata.WEC.Radius;
pd.WEC.Draft        = pdata.WEC.Draft;
pd.WEC.RD           = pdata.WEC.RD;
x                   = pdata.WEC.x;
y                   = pdata.WEC.y;
PTO_damping         = pdata.WEC.PTO_damping;
PTO_stiffness       = pdata.WEC.PTO_stiffness;

% Avoid plots
pd.General.plotflag = 0;

for i = 1:pdata.General.n_wec

    pd.WEC.PTO_damping   = PTO_damping(1,i);
    pd.WEC.PTO_stiffness = PTO_stiffness(1,i);
    pd.WEC.x             = x(i,1);
    pd.WEC.y             = y(i,1);
    pd                   = run_simulation(pd,'Q_factor');

    % Total expected power over 30 years for each device
    q_(i,1)              = sum(pd.Results.Power.P_total_vec); 
    
end

q_factor = sum(pdata.Results.Power.P_total_vec,2)./(q_);
Q_factor = sum(sum(pdata.Results.Power.P_total_vec,2))/sum(q_);

end