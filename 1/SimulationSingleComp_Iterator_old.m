clear all;
clc;
    
% Simulation parameters
p.question = 0;        %Question in assignment

p.t0 = 0;              %tspan start                        [ms]
p.tf = 25;             %%simulation time                   [ms]
p.h = .1;             %%time step
p.V0 = 0;              %%initial voltage                   [mV]

tspan=p.t0:p.h:p.tf;

%% Constants
p.C    = 1;            % membrane capacitance              [uF/cm^2]
p.gbar_K   = 36;       % max. potassium (K) conductance    [mS/cm^2]
p.gbar_Na  = 120;      % max. sodium (Na) conductance      [mS/cm^2]
p.gbar_L   = 0.3;      % leak (L) conductance              [mS/cm^2]

Ts=[50,22.2,6.3,0,-20];
cc=hsv(numel(Ts));
figure1=figure; 
hold on;

ylabel('Voltage (mV)')
xlabel('time (ms)')
title('How temprature effects Voltage over time in an artificial Neuron')
%axes1 = axes('Color',[0.313725501298904 0.313725501298904 0.313725501298904]);

for iT=1:numel(Ts)

    hold on
%% Calculating equibrilium (ion reversval potential) using nernst equation
p.T     = Ts(iT);
%p.T     = 22.2;        % temperature                       [Celsius]
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

%===plot Voltage===%


plot(tspan, V, 'color',cc(iT,:), 'LineWidth',1.5);

%axes([p.t0 p.tf -100 60]);
%{
set(gca,'XTick',p.t0:10:p.tf,...
        'YTick',-110:10:60);
%}


end



%hold(axes1,'all');
t = cellstr(num2str(Ts .'));
h=legend(t{:});
%set(h,'color','w');

BaseLine = ones(1,numel(tspan)) * p.RestVoltage;

hold off