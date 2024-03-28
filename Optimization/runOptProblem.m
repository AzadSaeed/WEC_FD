function [y,fval,pdata] = runOptProblem(pdata)

% This function solves the optimization problem based on the problem
% options, and the case study.
%
%
% Primary Contributor: Saeed Azad, PhD


switch upper(pdata.Opt.Solver)

    case 'FMINCON'

        [y,fval,pdata,history] = runfmincon(pdata);
        pdata.Results.Opt_history = history;
        close all;


    case 'SURROGATEOPT'

        [y,fval,pdata, exitflag, output,trials] = runsurrogateopt(pdata);
        pdata.Results.Opt_outout = output;
        pdata.Results.Opt_exitflag = exitflag;
        pdata.Results.Opt_trials = trials;

    case 'GA'

        [y,fval,pdata,exitflag,output,population,scores] = runga(pdata);
        pdata.Results.Opt_outout = output;
        pdata.Results.Opt_exitflag = exitflag;
        pdata.Results.Opt_populations = population;
        pdata.Results.Opt_scores = scores;
end



% Ensure the solution is a column
y = y(:);

end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [y,fval,pdata,history] = runfmincon(pdata)

% Set up shared variables with outfun
history.x             = [];
history.fval          = [];
history.firstorderopt = [];
history.funccount     = [];
history.fighandles    = [];


% Scaled, randomized starting point
if ~pdata.Opt.Hybrid
    if ~isfield(pdata.Opt,'rng')
        rng(5)
    else
        rng(pdata.Opt.rng)
    end
    x0_s                   = rand(pdata.Opt.nt,1);

    % Scale to obtain range and offset
    [~,pdata.Opt.D_scale,pdata.Opt.b_scale] = scale_var(x0_s,pdata.Opt.lb,pdata.Opt.ub);

elseif pdata.Opt.Hybrid

    x0_s = pdata.Opt.X_Init;
    [x0_s,pdata.Opt.D_scale,pdata.Opt.b_scale] = scale_var(x0_s,pdata.Opt.lb,pdata.Opt.ub);
    
end


% Scale to obtain range and offset
[~,pdata.Opt.D_scale,pdata.Opt.b_scale] = scale_var(x0_s,pdata.Opt.lb,pdata.Opt.ub);

% Define normalized bounds
LB                     = -1*ones(size(pdata.Opt.lb));
UB                     = ones(size(pdata.Opt.lb));

options = optimoptions('fmincon','Display',pdata.Opt.Display,'Algorithm',pdata.Opt.Algorithm,...
    'UseParallel',pdata.Opt.Parflag,'MaxIter',pdata.Opt.MaxFuncEval,...
    'OptimalityTolerance',pdata.Opt.Tolerance,'FiniteDifferenceStepSize',...
    pdata.Opt.FiniteDifferenceStepSize,...
    'FiniteDifferenceType','central','OutputFcn',@outfun1);

if pdata.Opt.SatFlag
    if ~pdata.General.OCFlag
        timerVal = tic;
        [y, fval] = fmincon(@(x)OBJ_Calculations(x,pdata),x0_s,[],[],[],...
            [],LB,UB,@(x)Constr(x,pdata),options);
        pdata.array_time = toc(timerVal);
    % else
    %     timerVal = tic;
    %     [y, fval] = Callfmincon(x0_s,LB, UB, options, pdata);
    %     pdata.array_time = toc(timerVal);
    end
else
    timerVal = tic;
    [y, fval] = fmincon(@(x)OBJ_Calculations(x,pdata),x0_s,[],[],[],[],...
        LB,UB,@(x)Constr(x,pdata),options);
    pdata.array_time = toc(timerVal);
end


function stop = outfun1(x,optimValues,state)
     stop = false;
 
     switch state

         case 'iter'

           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x'];
           history.firstorderopt = [history.firstorderopt;optimValues.firstorderopt];
           history.funccount     = [history.funccount; optimValues.funccount];

         case 'done'
             % fighandle1 = CreateOptHisPlots(history.x(:,i),'x_his');
             fighandle1 = CreateOptHisPlots(history,'fval_his');
             fighandle2 = CreateOptHisPlots(history,'funccount_hist');
             fighandle3 = CreateOptHisPlots(history,'firstorderopt_hist');
             history.fighandles = [ history.fighandles; fighandle1, fighandle2, fighandle3];

         otherwise
     end
end

Names = {'fval_his','funccount_hist', 'firstorderopt_hist'};
for i=1:length(history.fighandles)
    savename = strcat(pdata.General.plot_dir, filesep, strcat(Names{1,i},'.pdf'));
    exportgraphics(history.fighandles(1,i),savename);
end


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function [y,fval,pdata, exitflag, output,trials] = runsurrogateopt(pdata)


% Scaled, randomized starting point
if ~isfield(pdata.Opt,'rng')
    rng(5)
else
    rng(pdata.Opt.rng)
end
x0_s                   = rand(pdata.Opt.nt,1);

% Scale to obtain range and offset
[~,pdata.Opt.D_scale,pdata.Opt.b_scale] = scale_var(x0_s,pdata.Opt.lb,pdata.Opt.ub);

% Define normalized bounds
LB                     = -1*ones(size(pdata.Opt.lb));
UB                     = ones(size(pdata.Opt.lb));

% Initial condition
startpts               = x0_s;

if pdata.Opt.CheckPointFlag

    % Define checkpoint file name - expensive step
    Checkpointfilename = strcat('CheckP_',pdata.General.CaseStudy,'_',...
        string(pdata.General.n_wec),'WEC_',pdata.WaveData.Region,'.mat');

    if isfile(Checkpointfilename)
        timerVal = tic;
        [y,fval,exitflag,output,trials] = surrogateopt(Checkpointfilename);
        pdata.array_time = toc(timerVal);
    else

        % Surrogateopt options-plot: 'surrogateoptplot'
        options = optimoptions('surrogateopt','PlotFcn',pdata.Opt.PlotFunc,...
            'InitialPoints',startpts,'MaxFunctionEvaluations',pdata.Opt.MaxFuncEval, ...
            'UseParallel', pdata.Opt.Parflag,"CheckpointFile",Checkpointfilename);

        timerVal = tic;
        [y, fval, exitflag, output, trials] = surrogateopt(@(x)Objconstr_susurrogateopt(x,pdata),...
            LB,UB,[],[],[],[],[],options);
        pdata.array_time = toc(timerVal);

    end
else

    % Surrogateopt options
    options = optimoptions('surrogateopt','PlotFcn',pdata.Opt.PlotFunc,...
        'InitialPoints',startpts,'MaxFunctionEvaluations',pdata.Opt.MaxFuncEval, ...
        'UseParallel', pdata.Opt.Parflag);
    timerVal = tic;
    [y, fval, exitflag, output, trials] = surrogateopt(@(x)Objconstr_susurrogateopt(x,pdata),...
        LB,UB,[],[],[],[],[],options);
    pdata.array_time = toc(timerVal);
    
end



end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [y,fval,pdata,exitflag,output,population,scores] = runga(pdata)


if ~isfield(pdata.Opt,'rng')
    rng(5)
else
    rng(pdata.Opt.rng)
end

x0_s                   = rand(pdata.Opt.nt,1);

% Scale to obtain range and offset
[~,pdata.Opt.D_scale,pdata.Opt.b_scale] = scale_var(x0_s,pdata.Opt.lb,pdata.Opt.ub);

% Define normalized bounds
LB                     = -1*ones(size(pdata.Opt.lb));
UB                     = ones(size(pdata.Opt.lb));


% Define the Creation and mutation functions for GA
if strcmpi(pdata.General.CaseStudy,'C_OPT')

    Cre_Func = 'gacreationuniform';
    Mu_Func = 'mutationgaussian';
    Cro_Func = 'crossoverscattered';

    % Suppress warning for this case
    id = 'globaloptim:constrvalidate:unconstrainedMutationFcn';
    warning('off',id);

else
    Cre_Func = {'gacreationnonlinearfeasible','UseParallel',pdata.Opt.Parflag,'NumStartPts',pdata.Opt.PS};
    Mu_Func = 'mutationadaptfeasible';
    Cro_Func = {@crossoverintermediate, ones(1,pdata.Opt.nt)};
end





options = optimoptions('ga','PlotFcn',...
    {@gaplotbestf, @gaplotstopping},...
    'ConstraintTolerance', pdata.Opt.ConstraintTolerance, 'CrossoverFraction',...
    pdata.Opt.CrossoverFraction,'FunctionTolerance',pdata.Opt.FunctionTolerance,...
    'MaxGenerations',pdata.Opt.MG,...
    'MaxStallGenerations', pdata.Opt.MaxStallGen,'PopulationSize',pdata.Opt.PS, ...
    'UseParallel', pdata.Opt.Parflag,'CreationFcn',Cre_Func,...
    'SelectionFcn',pdata.Opt.Se_Func, 'MutationFcn',Mu_Func, 'CrossoverFcn',Cro_Func);

timerVal = tic;
[y,fval,exitflag,output,population,scores] = ga(@(x)Objconstr_GA(x,pdata),...
    pdata.Opt.nt,[],[],[],[],LB,UB,@(x)Constr(x,pdata),[],options);
pdata.array_time = toc(timerVal);

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function fighandle = CreateOptHisPlots(history,plotflag)

commonFigureProperties;
hf = figure;
hf.Color = 'w';
hold on;
iter  = 1:length(history.fval);
switch plotflag

    case 'x_his'
        
        plot(iter,history.x(:,1),'.-','LineWidth',linewidth, 'MarkerSize', 16, 'Color', C.blue(9,:))
        plot(iter,history.x(:,2),'.-','LineWidth',linewidth, 'MarkerSize', 16, 'Color', C.red(9,:))
        leg = {'\textrm{Normalized}~$B_{\textrm{pto}}$','$\textrm{Normalized}~k_{\textrm{pto}}$'};
        xL  = 'Iterations';
        yL  = '$\textrm{Optimization Variables}$';
        legend(leg,'Location','northoutside')  

    case 'fval_his'
        
        plot(iter,history.fval,'.-','LineWidth',linewidth, 'MarkerSize', 16, 'Color', C.purple(9,:))
        xL  = 'Iterations';
        yL  = '$\textrm{objective Function value}$';
        leg = {''};

    case 'funccount_hist'
        
        plot(iter,history.funccount,'.-','LineWidth',linewidth, 'MarkerSize', 16, 'Color', C.orange(9,:))
        xL  = 'Iterations';
        yL  = '$\textrm{Function Count}$';


    case 'firstorderopt_hist'

        plot(iter,history.firstorderopt,'.-','LineWidth',linewidth, 'MarkerSize', 16, 'Color', C.green(7,:))
        xL  = 'Iterations';
        yL  = '$\textrm{Function Count}$';

end

%title(titletext)
      
xlabel(xL, 'FontSize', fontlabel)
ylabel(yL, 'FontSize', fontlabel)

ax = gca;
ax.FontSize = fonttick;
taskflag = 'axes';
commonFigureTasks;
taskflag = 'legend';
commonFigureTasks;
fighandle = gca;

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
