function [result] = SimulationSingleCompartment(p)
%{
if nargin < 1
    p.question = 999;        %Question in assignment
end
%}

%% q5
if (p.question == 5) %TODO: put it in a place more suitable
    
    global lastNormRndCurrent;
    global nextNormRndCalcPercentile;

    lastNormRndCurrent=-inf;
    nextNormRndCalcPercentile=-inf;
end
    
    
% Simulation parameters
tspan=p.t0:p.h:p.tf;

%% Constants
p.V0 = 0;              %%initial voltage                   [mV]
p.C    = 1;            % membrane capacitance              [uF/cm^2]
p.gbar_K   = 36;       % max. potassium (K) conductance    [mS/cm^2]
p.gbar_Na  = 120;      % max. sodium (Na) conductance      [mS/cm^2]
p.gbar_L   = 0.3;      % leak (L) conductance              [mS/cm^2]

%% Calculating equibrilium (ion reversval potential) using nernst equation
p.T     = 22.2;        % temperature                       [Celsius]
%p.T    = 6.3;        % temperature                       [Celsius]
p.K_Out = 4;
p.K_In = 155 ;
p.K_Valance = 1;
p.Na_Out = 145;
p.Na_In = 12;
p.Na_Valance = 1;

%{
p.E_K   = -12;          % K equilibrium potential           [mV] 
p.E_Na  = 115;          % N equilibrium potential           [mV]
p.E_l   = 10.598;       % L equilibrium potential           [mV]
%}

%Parameters used maily for display
p.RestVoltage = -70;   % membrane resting potential        [mV]

p.E_K = nernst_fun(p.T, p.K_Valance, p.K_Out, p.K_In) ...
    -p.RestVoltage;
p.E_Na = nernst_fun(p.T, p.Na_Valance, p.Na_Out, p.Na_In) ...
    -p.RestVoltage;
p.E_l = -50.0 -p.RestVoltage;


%determines when to trace
global LastPercentileToTrace;
LastPercentileToTrace=-inf;

%% Initial state variable
n0 = nmh_Inf_fun(alpha_n_fun(p.V0), beta_n_fun(p.V0)); %Equation 9
m0 = nmh_Inf_fun(alpha_m_fun(p.V0), beta_m_fun(p.V0)); %Equation 9
h0 = nmh_Inf_fun(alpha_h_fun(p.V0), beta_h_fun(p.V0)); %Equation 9

[tspan,outArrSingle] = ode45( ...
        @(t,outArrSingle) hodgkinHuxleyDE(p,tspan,t, ...
                        outArrSingle), ...
        tspan, [p.V0 n0 m0 h0]);

%save('single.outArrSingle.mat','outArrSingle','-mat')

V = outArrSingle(:,1);
n = outArrSingle(:,2);
m = outArrSingle(:,3);
h = outArrSingle(:,4);

V = V+p.RestVoltage; %Set resting potential to -70mv intead of relative 0mv

result.p = p;
result.V = V;
result.m = m;
result.n = n;
result.h = h;
result.tspan = tspan;



