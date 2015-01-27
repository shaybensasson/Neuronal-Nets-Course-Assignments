clc;
clear all;
close all;

%% setup simulation parameters 
p.TRACE = 1;
p.PLOT_MULTI = 1;
p.RENDER_MOVIE = 1;
p.SAVE_MOVIE = 1;

p.tf     = 15;    % ms
p.dt    = .01; % ms
tspan  = 0:p.dt:p.tf;

%% axon properties
p.comps   = 20;                % number of compartments
%p.elec    = int8(p.comps/2);                % compartment index of stimulating electrode
p.elec    = 1;                % compartment index of stimulating electrode
p.RA      = 1e-6;               % specific intracellular resistivity
p.a       = 1e-3;             % compartment radius

p.L       = 300;      %axon length (cm)
p.l       = p.L./p.comps;          % compartment length (cm)
p.Ra      = (p.RA*p.l)/(2*pi*p.a.^2); % intracellular resistance


%% Stimulus
feedExternalCurrent = @(t, index) iif (t>=1 && t<=10, 10, 0);

%% HH Parameters
p.V0  = 0;      % mV
p.Cm      = 1;      % uF/cm2
p.gbar_Na = 120;    % mS/cm2
p.gbar_K  = 36;     % mS/cm2
p.gbar_l  = 0.3;    % mS/cm2

%% Calculating equibrilium (ion reversval potential) using nernst equation
p.T     = 22.2;        % temperature                       [Celsius]
% p.T    = 6.3;        % temperature                       [Celsius]
p.K_Out = 4;
p.K_In = 155 ;
p.K_Valance = 1;
p.Na_Out = 145;
p.Na_In = 12;
p.Na_Valance = 1;

%{
p.E_K   = -12;          % K equilibrium potential           [mV] 
p.E_Na  = 115;          % N equilibrium potential           [mV]
p.E_l   = 10.613;       % L equilibrium potential           [mV]
%}

% used mainly for display proposes
p.RestVoltage = -70;   % membrane resting potential        [mV]
p.PeakVoltage = +55;   % AP peak potential                 [mV]


p.E_K = nernst_fun(p.T, p.K_Valance, p.K_Out, p.K_In) ...
    -p.RestVoltage;
p.E_Na = nernst_fun(p.T, p.Na_Valance, p.Na_Out, p.Na_In) ...
    -p.RestVoltage;
p.E_l = -50.0 -p.RestVoltage;


%% Functions
% K channel
alpha_n = @(v) iif (v ~= 10, 0.01.*(-v + 10)./(exp((-v + 10)./10) - 1), 0.1);
%alpha_n = @(v) 0.01.*(-v + 10)./(exp((-v + 10)./10) - 1);
beta_n  = @(v) 0.125*exp(-v./80);
n_inf   = @(v) alpha_n(v)./(alpha_n(v) + beta_n(v));

% Na channel (activating)
%alpha_m = @(v) iif (v ~= 25, 0.1.*(-v + 25)./(exp((-v + 25)./10) - 1), 1);
alpha_m = @(v) 0.1.*(-v + 25)./(exp((-v + 25)./10) - 1);
beta_m  = @(v) 4*exp(-v./18);
m_inf   = @(v) alpha_m(v)./(alpha_m(v) + beta_m(v));

% Na channel (inactivating)
alpha_h = @(v) 0.07*exp(-v./20);
beta_h  = @(v) 1./(exp((-v + 30)./10) + 1);
h_inf   = @(v) alpha_h(v)./(alpha_h(v) + beta_h(v));

%% allocate state variables
VmAll      = zeros(p.comps,numel(tspan));  % mV
VmAll(:,1) = p.V0;

mAll       = ones(p.comps,1)*m_inf(p.V0);
hAll       = ones(p.comps,1)*h_inf(p.V0);
nAll       = ones(p.comps,1)*n_inf(p.V0);


%% connection matrix
conMat = zeros(p.comps,p.comps);
for comp=1:p.comps
  if comp == 1 %first row - end compartment - right voltage drop only
      conMat(comp,1:2)     = [1 -1];
  elseif comp == p.comps %last row - end compartment
      conMat(comp,p.comps-1:p.comps)   = [-1 1]; %left voltage drop only
  else
      conMat(comp,comp-1:comp+1) = [-1 2 -1]; %compartment ion flow in both ways
  end
end

%% simulate model
[tspan,outArr] = ode45( ...
        @(t,outArr) hodgkinHuxleyMultiCompDE(p,tspan,t, ...
                        conMat, outArr, ...
                        feedExternalCurrent), ...
        tspan, [VmAll(:,1);nAll;mAll;hAll]);

outArr = reshape(outArr, [numel(tspan) p.comps 4]);

%save('multi.outArr.mat','outArr','-mat')

%seperate param from other params
VmAll = outArr(:,:,1); 
nAll = outArr(:,:,2);
mAll = outArr(:,:,3);
hAll = outArr(:,:,4);

VmAll = VmAll+p.RestVoltage; %Set resting potential to -70mv intead of relative 0mv


%% Calculating Velocity
%====================================
start_d = 1;
start_t = 1;
for comp=1:1:p.comps
    [m1,i1] = max(VmAll(:,comp));
    if (m1>0)
        start_d = comp;
        start_t = i1;
        break
    end
end


end_d = 1;
end_t = 1;
for comp=p.comps:-1:1
    [m2,i2] = max(VmAll(:,comp));
    if (m2>0)
        end_d = comp;
        end_t = i2;
        break
    end
end

%take the length and divide it by time elapsed and normalize to units of ms
range_l = abs(end_d-start_d)*p.l;
range_l = range_l .* 100; %convert cm to m
range_t = abs(end_t-start_t)*p.dt .* 1000; %convert ms to s
velocity = 10*(range_l ./ range_t);

%{
fprintf('p.RA: %s\n', num2str(p.RA,'%1.2e'));
fprintf('p.a: %s cm\n', num2str(p.a,'%1.2e'));
fprintf('p.L: %s cm\n', num2str(p.L,'%1.2e'));
%}

fprintf('VELOCITY: %s m/s\n', num2str(velocity,10));

    
%% plot membrane potential by compartment
countFigures=0;
if (p.PLOT_MULTI)
    SUBPLOTS_PER_FIGURE = 10;
    bulks=ceil(p.comps/SUBPLOTS_PER_FIGURE);
    for comp=1:p.comps
        bulkCompIndex = mod(comp,SUBPLOTS_PER_FIGURE);
        bulkCompIndex = iif(bulkCompIndex==0, ...
            SUBPLOTS_PER_FIGURE,bulkCompIndex);
        
        if (bulkCompIndex == 1)
            countFigures = countFigures+1;
            fig = figure(countFigures);
        end
        
        subplot(SUBPLOTS_PER_FIGURE,1,bulkCompIndex);
        plot(tspan, VmAll(:,comp),'linewidth',1.5);
        title(sprintf('Compartment #%d', comp));
        if (bulkCompIndex==1)
          ylabel('Membrane Potential (mV)');
        end

        if (bulkCompIndex==SUBPLOTS_PER_FIGURE)
          xlabel('Time (msec)');
        end

        axis([0 p.tf -110 60]);
    end
    
end
    
    
    
    %saveas(fig,'MultiCompAPTravel','jpg');
%else
    if (p.RENDER_MOVIE)
        countFigures = countFigures+1;
        fig = figure(countFigures);
        rect = get(fig,'Position'); 
        rect(1:2) = [0 0]; 

        set(fig,'DoubleBuffer','on');
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');

        if (p.SAVE_MOVIE)
            writerObj = VideoWriter('MultiCompAPTravel.avi');
            writerObj.FrameRate=5; %5 fps
            open(writerObj);
        end



        for comp=1:p.comps
          h = plot(tspan, VmAll(:,comp),'linewidth',1.5);
          %movegui(h, 'onscreen');
          title(sprintf('Compartment #%d', comp));
          ylabel('Membrane Potential (mV)');
          xlabel('Time (msec)');
          axis([0 p.tf -110 60]);
          text(1,1,['VELOCITY:' num2str(velocity,10) 'm/s'],...
            'HorizontalAlignment','left');

          drawnow; 
          if (p.SAVE_MOVIE) 
              writeVideo(writerObj,getframe(gcf,rect));
          end
          pause(.1);
        end

        if (p.SAVE_MOVIE)
            for i=1:15 %write another 3 seconds suffix
                writeVideo(writerObj,getframe(gcf,rect));
            end
            
            close(writerObj);
        end
       
    end
    
%end

%VELOCITY: 20.51048314 m/s
