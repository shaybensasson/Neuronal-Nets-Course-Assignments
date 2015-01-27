tic;
clc;
clear all;
close all;

%% simulation paramters
p.t0 = 0;              %tspan start                        [ms]
p.tf = 1000;             %%simulation time                   [ms]
p.h = 1;             %%time step


MIN_CURRENT = 0;
STEP_SIZE=1;
MAX_CURRENT = 80;

PLOT_VOLTAGE = 0;


internalCounter=0;

for externalCurrent=MIN_CURRENT:STEP_SIZE:MAX_CURRENT
    internalCounter=internalCounter+1;
    p.question = 999;        %Question in assignment
    p.NoQuestionConstantCurrent = externalCurrent;
    
    %{
    p.question = 5;
    p.Mu = externalCurrent;
    %p.Sigma = 7;
    p.Sigma = 10;
    %}
        
    result = SimulationSingleCompartment(p);
    
    pks = findpeaks(result.V,'MINPEAKHEIGHT',0);
    %CurrentLevel
    Acc(internalCounter,1) = externalCurrent;
    %APNumber
    countPks = size(pks,1);
    Acc(internalCounter,2) = countPks;
    
    fprintf('externalCurrent(%dmA): APNumber %d\n', externalCurrent, ...
        countPks);
    
    
    %===plot Voltage===%
    if (PLOT_VOLTAGE)
        hFig=figure(internalCounter);
        subplot(2,1,1);
        hold on
        plot(result.tspan,result.V,'LineWidth',2);
        BaseLine = ones(1,numel(result.tspan)) * result.p.RestVoltage;
        plot(result.tspan,BaseLine,':k', 'LineWidth',1)
        ylabel('Voltage (mV)')
        xlabel('time (ms)')
        if (p.question==999)
            title(sprintf('Voltage over time in an artificial Neuron, Constant I=%d', ...
                externalCurrent));
        elseif (p.question==5)
            title(sprintf('Mu:%d, Sigma:%d, tf:%d ms', p.Mu,p.Sigma,p.tf));
        else
            title('Voltage over time in an artificial Neuron')
        end
        axis([p.t0 p.tf -100 60]);


        set(gca,'XTick',p.t0:floor(p.tf/10):p.tf,...
                'YTick',-110:10:60);


        hold off

        subplot(2,1,2);
        plot(result.tspan,result.n, ':b');
        hold on
        plot(result.tspan,result.m, 'r');
        hold on
        plot(result.tspan,result.h, 'g--');
        legend('n gates (K+)','m gates (Na+)', 'h gates (Na+)');
        xlabel('time');
        title('Gates over time');
        hold off

        %saveas(hFig, sprintf('Mu%d_Sigma%d_t%d.png', p.Mu,p.Sigma,p.tf));
    end

    
end

%IAA(x,1) = x;
%IAA(x,2) = AANumber;




t=toc;
disp(['Elapsed time is ' datestr(datenum(0,0,0,0,0,t),'HH:MM:SS') '.']);

beep;

