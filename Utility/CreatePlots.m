function CreatePlots(pdata)

switch upper(pdata.General.CaseStudy)
    
    case 'P_OPT'
        Create_P_opt_plots(pdata);
    case 'C_OPT'
        Create_C_opt_plots(pdata);
    case 'L_OPT'
        Create_L_opt_plots(pdata);
    case 'LP_OPT'
        Create_LP_opt_plots(pdata)
    case 'PC_OPT'
        Create_PC_opt_plots(pdata)
    case 'PCL_OPT'
        Create_PCL_opt_plots(pdata)
end
end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Create_P_opt_plots(pdata)

% Plot the Prescribed layout
commonFigureProperties;
hf = figure;
hf.Color = 'w';
plot(pdata.Results.x, pdata.Results.y,'o','MarkerFaceColor', C.blue(7,:),...
    'MarkerEdgeColor', C.yellow(7,:),'markerSize', 8)
xlabel('$\textrm{x}~[m]$')
ylabel('$\textrm{y}~[m]$')
title('Prescribed layout')
plot_name = strcat(pdata.General.Location_Flag,'_Prescribed_layout','.pdf');
exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
close(hf)

% If irregular, plot power generated in each year
if strcmpi(pdata.General.WaveType,'IRREGULAR')

    Volume = pi*pdata.Results.Radius^2*pdata.Results.Draft*2;
    
    % Power in kW
    Power  = sum(pdata.Results.Power.P_total_vec)/1000;
    PV     = Power./Volume;

    commonFigureProperties;
    hf = figure;
    hf.Color = 'w';
    
    plot(Power,'o-','MarkerFaceColor', C.blue(7,:),...
        'MarkerEdgeColor', C.blue(7,:),'markerSize', 6,'linewidth',linewidth)
    xlabel('$\textrm{Year}$')
    ylabel('$\textrm{Power}~[kW]$')

    yyaxis right
    plot(PV,'o-','MarkerFaceColor', C.green(7,:),...
        'MarkerEdgeColor', C.green(7,:),'markerSize', 6,'linewidth',linewidth)
    ylabel('$\textrm{Power per Volume}~[kW/m^3]$')

    legend({'Power','Power per Volume'},'Location','best')
    plot_name = strcat('Annual_Power','.pdf');
    exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
    close(hf)
end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Create_C_opt_plots(pdata)

% Plot the Prescribed layout
commonFigureProperties;
hf = figure;
hf.Color = 'w';
plot(pdata.Results.x, pdata.Results.y,'o','MarkerFaceColor', C.blue(7,:),...
    'MarkerEdgeColor', C.yellow(7,:),'markerSize', 8)
xlabel('$\textrm{x}~[m]$')
ylabel('$\textrm{y}~[m]$')
title('Prescribed layout')
plot_name = strcat(pdata.General.Location_Flag,'_Prescribed_layout','.pdf');
exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
close(hf)

% If irregular, plot power generated in each year
if strcmpi(pdata.General.WaveType,'IRREGULAR')

    Volume = pi*pdata.Results.Radius^2*pdata.Results.Draft*2;
    
    % Power in kW
    Power  = sum(pdata.Results.Power.P_total_vec)/1000;
    PV     = Power./Volume;

    commonFigureProperties;
    hf = figure;
    hf.Color = 'w';
    
    plot(Power,'o-','MarkerFaceColor', C.blue(7,:),...
        'MarkerEdgeColor', C.blue(7,:),'markerSize', 6,'linewidth',linewidth)
    xlabel('$\textrm{Year}$')
    ylabel('$\textrm{Power}~[kW]$')

    yyaxis right
    plot(PV,'o-','MarkerFaceColor', C.green(7,:),...
        'MarkerEdgeColor', C.green(7,:),'markerSize', 6,'linewidth',linewidth)
    ylabel('$\textrm{Power per Volume}~[kW/m^3]$')

    legend({'Power','Power per Volume'},'Location','best')
    plot_name = strcat('Annual_Power','.pdf');
    exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
    close(hf)
end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Create_L_opt_plots(pdata)

% Safety sistance constraints
SF      = pdata.General.SafeDist;
centers = [pdata.Results.x,pdata.Results.y];
radii   = ((2*pdata.Results.Radius) + SF)/2;

% Farm area
farm_area = 0.5*sqrt(20000*pdata.General.n_wec);
FA_x        = [0,0,farm_area,farm_area, 0];
FA_y        = [-farm_area,farm_area,farm_area, -farm_area, -farm_area];

% Plot the optimized layout
commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on
plot(pdata.Results.x, pdata.Results.y,'o','MarkerFaceColor', C.blue(7,:),...
    'MarkerEdgeColor', C.yellow(7,:),'markerSize', 8)

viscircles(centers,repmat(radii,[length(centers),1]),'Color', C.red(9,:),...
    'LineStyle',':','LineWidth',2)

plot(FA_x,FA_y,':','Color',C.red(5,:),'markerEdgeColor',C.red(5,:),'markerSize',3)
axis equal


% Label the WECs
for i = 1:pdata.General.n_wec
       wec_id = string(i);
       text(pdata.Results.x(i,1), pdata.Results.y(i,1)+8,wec_id);
end

xlabel('$\textrm{x}~[m]$')
ylabel('$\textrm{y}~[m]$')
title('Optimized layout')
plot_name = strcat('optimized_layout','.pdf');
exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
close(hf)

% If irregular, plot power generated in each year
if strcmpi(pdata.General.WaveType,'IRREGULAR')

    Volume = pi*pdata.Results.Radius^2*pdata.Results.Draft*2;

    % Power in kW
    Power  = sum(pdata.Results.Power.P_total_vec)/1000;
    PV     = Power./Volume;

    commonFigureProperties;
    hf = figure;
    hf.Color = 'w';
    
    plot(Power,'o-','MarkerFaceColor', C.blue(7,:),...
        'MarkerEdgeColor', C.blue(7,:),'markerSize', 6,'linewidth',linewidth)
    xlabel('$\textrm{Year}$')
    ylabel('$\textrm{Power}~[kW]$')

    yyaxis right
    plot(PV,'o-','MarkerFaceColor', C.green(7,:),...
        'MarkerEdgeColor', C.green(7,:),'markerSize', 6,'linewidth',linewidth)
    ylabel('$\textrm{Power per Volume}~[kW/m^3]$')

    legend({'Power','Power per Volume'},'Location','best')
    plot_name = strcat('Annual_Power','.pdf');
    exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
    close(hf)
end

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Create_LP_opt_plots(pdata)

% Safety sistance constraints
SF      = pdata.General.SafeDist;
centers = [pdata.Results.x,pdata.Results.y];
radii   = ((2*pdata.Results.Radius) + SF)/2;

% Farm area
farm_area = 0.5*sqrt(20000*pdata.General.n_wec);
FA_x        = [0,0,farm_area,farm_area, 0];
FA_y        = [-farm_area,farm_area,farm_area, -farm_area, -farm_area];

% Plot the optimized layout
commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on
plot(pdata.Results.x, pdata.Results.y,'o','MarkerFaceColor', C.blue(7,:),...
    'MarkerEdgeColor', C.yellow(7,:),'markerSize', 8)

viscircles(centers,repmat(radii,[length(centers),1]),'Color', C.red(9,:),...
    'LineStyle',':','LineWidth',2)

plot(FA_x,FA_y,':','Color',C.red(5,:),'markerEdgeColor',C.red(5,:),'markerSize',3)
axis equal


% Label the WECs
for i = 1:pdata.General.n_wec
       wec_id = string(i);
       text(pdata.Results.x(i,1), pdata.Results.y(i,1)+8,wec_id);
end

xlabel('$\textrm{x}~[m]$')
ylabel('$\textrm{y}~[m]$')
title('Optimized layout')
plot_name = strcat('optimized_layout','.pdf');
exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
close(hf)

% If irregular, plot power generated in each year
if strcmpi(pdata.General.WaveType,'IRREGULAR')

    Volume = pi*pdata.Results.Radius^2*pdata.Results.Draft*2;

    % Power in kW
    Power  = sum(pdata.Results.Power.P_total_vec)/1000;
    PV     = Power./Volume;

    commonFigureProperties;
    hf = figure;
    hf.Color = 'w';
    
    plot(Power,'o-','MarkerFaceColor', C.blue(7,:),...
        'MarkerEdgeColor', C.blue(7,:),'markerSize', 6,'linewidth',linewidth)
    xlabel('$\textrm{Year}$')
    ylabel('$\textrm{Power}~[kW]$')

    yyaxis right
    plot(PV,'o-','MarkerFaceColor', C.green(7,:),...
        'MarkerEdgeColor', C.green(7,:),'markerSize', 6,'linewidth',linewidth)
    ylabel('$\textrm{Power per Volume}~[kW/m^3]$')

    legend({'Power','Power per Volume'},'Location','best')
    plot_name = strcat('Annual_Power','.pdf');
    exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
    close(hf)
end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Create_PC_opt_plots(pdata)

% Plot the Prescribed layout
commonFigureProperties;
hf = figure;
hf.Color = 'w';
plot(pdata.Results.x, pdata.Results.y,'o','MarkerFaceColor', C.blue(7,:),...
    'MarkerEdgeColor', C.yellow(7,:),'markerSize', 8)
xlabel('$\textrm{x}~[m]$')
ylabel('$\textrm{y}~[m]$')
title('Prescribed layout')
plot_name = strcat(pdata.General.Location_Flag,'_Prescribed_layout','.pdf');
exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
close(hf)

% If irregular, plot power generated in each year
if strcmpi(pdata.General.WaveType,'IRREGULAR')

    Volume = pi*pdata.Results.Radius^2*pdata.Results.Draft*2;
    
    % Power in kW
    Power  = sum(pdata.Results.Power.P_total_vec)/1000;
    PV     = Power./Volume;

    commonFigureProperties;
    hf = figure;
    hf.Color = 'w';
    
    plot(Power,'o-','MarkerFaceColor', C.blue(7,:),...
        'MarkerEdgeColor', C.blue(7,:),'markerSize', 6,'linewidth',linewidth)
    xlabel('$\textrm{Year}$')
    ylabel('$\textrm{Power}~[kW]$')

    yyaxis right
    plot(PV,'o-','MarkerFaceColor', C.green(7,:),...
        'MarkerEdgeColor', C.green(7,:),'markerSize', 6,'linewidth',linewidth)
    ylabel('$\textrm{Power per Volume}~[kW/m^3]$')

    legend({'Power','Power per Volume'},'Location','best')
    plot_name = strcat('Annual_Power','.pdf');
    exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
    close(hf)
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function Create_PCL_opt_plots(pdata)

% Safety sistance constraints
SF      = pdata.General.SafeDist;
centers = [pdata.Results.x,pdata.Results.y];
radii   = ((2*pdata.Results.Radius) + SF)/2;

% Farm area
farm_area = 0.5*sqrt(20000*pdata.General.n_wec);
FA_x        = [0,0,farm_area,farm_area, 0];
FA_y        = [-farm_area,farm_area,farm_area, -farm_area, -farm_area];

% Plot the optimized layout
commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on
plot(pdata.Results.x, pdata.Results.y,'o','MarkerFaceColor', C.blue(7,:),...
    'MarkerEdgeColor', C.yellow(7,:),'markerSize', 8)

viscircles(centers,repmat(radii,[length(centers),1]),'Color', C.red(9,:),...
    'LineStyle',':','LineWidth',2)

plot(FA_x,FA_y,':','Color',C.red(5,:),'markerEdgeColor',C.red(5,:),'markerSize',3)
axis equal


% Label the WECs
for i = 1:pdata.General.n_wec
       wec_id = string(i);
       text(pdata.Results.x(i,1), pdata.Results.y(i,1)+8,wec_id);
end

xlabel('$\textrm{x}~[m]$')
ylabel('$\textrm{y}~[m]$')
title('Optimized layout')
plot_name = strcat('optimized_layout','.pdf');
exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
close(hf)

% If irregular, plot power generated in each year
if strcmpi(pdata.General.WaveType,'IRREGULAR')

    Volume = pi*pdata.Results.Radius^2*pdata.Results.Draft*2;

    % Power in kW
    Power  = sum(pdata.Results.Power.P_total_vec)/1000;
    PV     = Power./Volume;

    commonFigureProperties;
    hf = figure;
    hf.Color = 'w';
    
    plot(Power,'o-','MarkerFaceColor', C.blue(7,:),...
        'MarkerEdgeColor', C.blue(7,:),'markerSize', 6,'linewidth',linewidth)
    xlabel('$\textrm{Year}$')
    ylabel('$\textrm{Power}~[kW]$')

    yyaxis right
    plot(PV,'o-','MarkerFaceColor', C.green(7,:),...
        'MarkerEdgeColor', C.green(7,:),'markerSize', 6,'linewidth',linewidth)
    ylabel('$\textrm{Power per Volume}~[kW/m^3]$')

    legend({'Power','Power per Volume'},'Location','best')
    plot_name = strcat('Annual_Power','.pdf');
    exportgraphics(gca,strcat(pdata.General.plot_dir,filesep,plot_name));
    close(hf)
end

end


