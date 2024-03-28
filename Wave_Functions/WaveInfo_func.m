function pdata = WaveInfo_func(pdata)

Wave_dir = strcat(pdata.General.pathstr,filesep, 'Wave_Functions');

switch pdata.WaveData.WDP

    case 'historical'

        FileName_general = strcat(Wave_dir,filesep,pdata.WaveData.Region,filesep,'USGS_Mat', filesep, 'USGS_',pdata.WaveData.site,'_gfdlrcp_hist');
        h_hist = load(strcat(FileName_general,'_height.mat'));
        t_hist = load(strcat(FileName_general,'_period.mat'));
        years = 30;

        series_h = (h_hist.height);
        series_t = (t_hist.period);

    case 'mid'

        FileName_general = strcat(Wave_dir,filesep,pdata.WaveData.Region,filesep,'USGS_Mat', filesep, 'USGS_',pdata.WaveData.site,'_gfdlrcp',pdata.WaveData.climate,'_mid');
        h_hist = load(strcat(FileName_general,'_height.mat'));
        t_hist = load(strcat(FileName_general,'_period.mat'));
        years = 12;


    case 'old_data'

        h_hist = load(strcat(Wave_dir,filesep,'Ex',filesep,'Mat_Files',filesep,'USGS_PI14_gfdlrcp85_historical','.mat'));
        t_hist = load(strcat(Wave_dir,filesep,'Ex',filesep,'Mat_Files',filesep,'USGS_PI14_gfdlrcp85_historical_mperiod','.mat'));
        years = 30;
        series_h = (h_hist.var);
        series_t = (t_hist.var);

end


t_len = years;
bin = 24*365;

h_mean = zeros(t_len,1);
t_mean = zeros(t_len,1);
h_var = zeros(t_len,1);
t_var = zeros(t_len,1);
cov_ht = {t_len,1};
cov_off = zeros(t_len,1);
h_yearly_pts = {t_len,1};
t_yearly_pts = {t_len,1};

for i = 1:t_len

    st = 1 + (i-1)*bin;
    fn = (i*bin);
    h_mean(i) = mean(series_h(st:1:fn));
    t_mean(i) = mean(series_t(st:1:fn));
    h_var(i) = var(series_h(st:1:fn));
    t_var(i) = var(series_t(st:1:fn));
    cov_ht{i} = cov(series_h(st:1:fn),series_t(st:1:fn));
    cov_off(i) = cov_ht{i}(1,2);
    h_yearly_pts{i} = series_h(st:1:fn);
    t_yearly_pts{i} = series_t(st:1:fn);

end


%% Define limits of integration

num_years = t_len;

max_t = zeros(num_years,1);
max_h = zeros(num_years,1);
min_t = zeros(num_years,1);
min_h = zeros(num_years,1);

for y = 1:num_years

    max_t(y) = max(t_yearly_pts{1,y});
    min_t(y) = min(t_yearly_pts{1,y});

    max_h(y) = max(h_yearly_pts{1,y});
    min_h(y) = min(h_yearly_pts{1,y});

end


% Defined domain of integration in y
% Using min(min) and max(max) still cuts the distribution, so adjust by 0.5
pdata.WaveData.a = min(min_h) - 0.5;
pdata.WaveData.b = max(max_h) + 0.5;


% Define domain of integration in x
pdata.WaveData.c = min(min_t) - 0.5;
pdata.WaveData.d = max(max_t) + 0.5;
 

%% Gauss quadrature 
pdata.WaveData.num_years = num_years;

% Create a surrogate model for probability distribution
pdata.WaveData  = gauss_q_pts_new(pdata.WaveData);


p_mat           = joint_prob_dist(t_yearly_pts,h_yearly_pts,pdata.WaveData,pdata);
pdata.WaveData.Pr_mat  = p_mat;

end