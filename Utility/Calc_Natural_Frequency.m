function x = Calc_Natural_Frequency(varargin)

if nargin ==1

    pdata  = varargin{1,1};
    K_pto  = pdata.Results.PTO_stiffness;
    Radius = pdata.Results.Radius;
    Draft  = pdata.Results.Draft;
    S      = pi*(Radius)^2;
    G      = pdata.General.rho*pdata.General.g*S;
    Mass   = pdata.General.rho*S*Draft*2;
    In     = pdata;


else

    In      = MS_input;
    Radius  = 0.5;
    Draft   = 0.5/1;
    S       = pi*(Radius)^2;
    G       = In.General.rho*In.General.g*S;
    Mass    = In.General.rho*S*Draft*2;
    In.General.n_wec = 1;
    In.Results.x = 0;
    In.Results.y = 0;
    K_pto   = 80000;
    pdata.Opt.CtrlFlag = 'FARM';

end


if In.General.n_wec==1
    x = SOlveEQ(Radius,Draft,Mass,K_pto,G,In);
else

    x = NaN(1,pdata.General.n_wec);
    
    for i=1:pdata.General.n_wec
        x(1,i) = SOlveEQ(Radius,Draft,Mass,K_pto(1,i),G,In, i);

        % Visualization
        % wn = linspace(0.1, 10,80);
        % for j=1:length(wn)
        %     P1(j) = (wn(j)^2)*(Mass + Calc_AddedMass(wn(j),Radius,Draft,In,i));
        % end
        % P2 = (K_pto(1,i)+G);
        % 
        % plot(wn,P1-P2, '.-r')
        % hold on
        % plot(wn,zeros(size(wn)),'.-b')
        % close all

    end

end

end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function x = SOlveEQ(Radius,Draft,Mass,K_pto,G,In, varargin)

if ~isempty(varargin)
    i = varargin{1,1}; % counter
else
    i = 1;
end



Method = 'fzero';

% Options
options = optimset('TolFun',10^-8,'MaxIter', 400,'Display','iter');

% Equation to solve
fun = @(wn) (wn^2)*(Mass + Calc_AddedMass(wn,Radius,Draft,In,i)) - (K_pto+G);

switch Method

    case 'fsolve'


        % starting point
        x0 = rand(1,1);
   

        % Solve the equation
        [x,~,exitflag,~]  = fsolve(fun,x0,options);

        if exitflag>0
            formatspec             = 'Equation solved. Natural Frequency is = %5.6f\n';
            fprintf(formatspec,x)
        elseif exitflag ==0
            formatspec             = 'Equation not solved. Number of iterations exceeded.\n';
            fprintf(formatspec)
        else
            formatspec             = 'Equation not solved. Number of iterations exceeded.\n';
            fprintf(formatspec)
        end

    case 'fzero'

        % starting point
        x0 = [10^-3 10];

        % Solve the equation
        [x,~,exitflag,~]  = fzero(fun,x0,options);

        if exitflag>0
            formatspec             = 'Equation solved. Natural Frequency is = %5.6f\n';
            fprintf(formatspec,x)
        elseif exitflag <0
            formatspec             = 'Equation not solved. Number of iterations exceeded.\n';
            fprintf(formatspec)
        end
end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function A = Calc_AddedMass(wn,Radius,Draft,In,i)

n_wec     = In.General.n_wec;
x         = In.Results.x;
y         = In.Results.y;
depth     = In.General.depth;
g         = In.General.g;
rho       = In.General.rho;
In.WEC.w  =  wn;

A_matrix = MyMSCoefficients(n_wec,[x,y],length(wn),Radius,Draft, depth, g, rho, wn,In.MS);

A = sum(A_matrix,2);
A = A(i,1);

% K_pto  = In.Results.PTO_stiffness;
% Radius = In.Results.Radius;
% Draft  = In.Results.Draft;
% S      = pi*(Radius)^2;
% G      = In.General.rho*In.General.g*S;
% Mass   = In.General.rho*S*Draft*2;
% val    = (wn^2)*(Mass + A) - (K_pto+G);
% 
% plot(wn,val,'.-r')
% hold on


end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function In = MS_input

% MS settings
In.MS.maxs  = 5;                          % order of interaction
In.MS.maxm  = 5;                          % order of series truncation for main fluid
In.MS.maxn  = 40;                         % order of series truncation for fluid below body
In.MS.maxj  = 15;                         % order of series truncation for diffraction problem
In.MS.fmaxw = 2;                          % max frequency in rad/s (for wave.m)
In.MS.fresw = 1;                          % number of frequencies for wave.m evaluations (linearly spaced from 0 to fmaxw)
In.MS.fres  = 50;

% Water characteristics
In.General.depth = 50;             % water depth (m)
In.General.g     = 9.81;           % gravity (m/s^2)
In.General.rho   = 1025;           % water density (kg/m^3)


end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
