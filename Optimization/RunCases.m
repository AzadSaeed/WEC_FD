function pdata = RunCases(pdata)
% This function is created to run all of the possible optimization studies.
% The cases are classified based on the level of complexity, and involve
% plant optimization 'P_opt', control optimization 'C_opt', layout
% optimization 'L_opt', plant and control optimization 'PC_opt', layout and
% plant optimization 'LP_opt', and plant, control, and layout optimization
% 'PCL_opt'
%
%
%
% Primary Contributer: Saeed Azad, PhD


switch upper(pdata.General.CaseStudy)

    case 'PCL_OPT'

        pdata = PCL_OPT_function(pdata);

    case 'C_OPT'

        pdata = C_OPT_function(pdata);

    case 'L_OPT'

        pdata = L_OPT_function(pdata);

    case 'P_OPT'

        pdata = P_OPT_function(pdata);

    case 'LP_OPT'

        pdata = LP_OPT_function(pdata);

    case 'PC_OPT'

        pdata = PC_OPT_function(pdata);


end


CreatePlots(pdata);

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = P_OPT_function(pdata)

Opt_Time               = tic;
[y,fval,pdata]         = runOptProblem(pdata);
pdata.Results.Opt_Time = toc(Opt_Time);
formatspec             = 'Plant optimization computational time is = %3.2f seconds\n';
fprintf(formatspec,pdata.array_time)
pdata = PostProcess(y,fval,pdata);
close all;

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = C_OPT_function(pdata)

Opt_Time               = tic;
[y,fval,pdata]         = runOptProblem(pdata);
pdata.Results.Opt_Time = toc(Opt_Time);
formatspec             = 'Control optimization computational time is = %3.2f seconds\n';
fprintf(formatspec,pdata.array_time)
pdata = PostProcess(y,fval,pdata);

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = L_OPT_function(pdata)


Opt_Time               = tic;
[y,fval,pdata]         = runOptProblem(pdata);
pdata.Results.Opt_Time = toc(Opt_Time);
formatspec             = 'Layout optimization computational time is = %3.2f seconds\n';
fprintf(formatspec,pdata.array_time)
pdata = PostProcess(y,fval,pdata);

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = LP_OPT_function(pdata)

Opt_Time              = tic;
[y,fval,pdata]         = runOptProblem(pdata);
pdata.Results.Opt_Time = toc(Opt_Time);
formatspec             = 'Layout and Plant optimization computational time is = %3.2f seconds\n';
fprintf(formatspec,pdata.array_time)
pdata = PostProcess(y,fval,pdata);


end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = PC_OPT_function(pdata)

Opt_Time              = tic;
[y,fval,pdata]         = runOptProblem(pdata);
pdata.Results.Opt_Time = toc(Opt_Time);
formatspec             = 'Plant and control optimization computational time is = %3.2f seconds\n';
fprintf(formatspec,pdata.array_time)
pdata = PostProcess(y,fval,pdata);


end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = PCL_OPT_function(pdata)

Opt_Time                    = tic;
[y,fval,pdata]              = runOptProblem(pdata);
pdata.Results.Opt_Time      = toc(Opt_Time);
formatspec                  = 'Plant, control, and layout optimization computational time is = %3.2f seconds\n';
fprintf(formatspec,pdata.array_time)
pdata = PostProcess(y,fval,pdata);

end

