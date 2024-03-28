function pdata = ProblemOptions(pdata,varargin)
% This functions enables the assignment and definition of various setting,
% options, and problem parameters.
%
% The varargin argument, if not empty, can have the following parameters in
% order: Region, Hydrodynamic Flag, Controller Flag, Initial optimization
%        RNG, Prescribed Layout, Solver, Power Saturation Limit, WaveType
% 
%
%
% The general options include: water density, gravitational constant, water 
% depth, safety distance for farms, number of frequencies, minimum and 
% maximum frequencies, plot flag, power conversion chain efficiency, operational
% availability, transmission efficiency, number of years, and Saturation 
% Method, which is set to PowerBound in this release. The wave type is also
% assigned in the general options 
%
%
% Optimization_options defines optimization-specific parameters, including
% bounds, number of variables
%
%
% Primary Contributer: Saeed Azad, PhD



if length(varargin) == 7
    pdata = General_options(pdata,varargin{1,7});
elseif length(varargin)==8
    pdata = General_options(pdata,varargin{1,7},varargin{1,8});
else
    pdata = General_options(pdata);
end


% Optimization options
if ~strcmpi(pdata.General.CaseStudy,'Sim')
    pdata = Optimization_options(pdata,varargin);
end

% Assign plot and solution directorie based on the study
if ~strcmpi(pdata.General.CaseStudy,'SIM')
    pdata = dir_(pdata);
end

% Wave Options
pdata = Wave_options(pdata,varargin{1,1});

% Multi scattering options
pdata = MS_options(pdata);


end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%








