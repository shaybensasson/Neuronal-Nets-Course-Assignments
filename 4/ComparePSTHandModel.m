close all;

ConstantsHeader();

%on question 3, we where asked to compare PSTH to model in Rep only
MODE = 'NonRep';
if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)) || ...
        Simulation.Phase < CONSTANTS.PHASES.GENERATOR)
    clearvars -except MODE;
    load(['AfterGenerator_' MODE '.mat'])
    ConstantsHeader();
end

NEURONS = length(Simulation.Neuron);
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;
SECONDS_OF_RATE_TO_DISPLAY = Simulation.SECONDS_OF_RATE_TO_DISPLAY;

%% plotting stims and after filters stims against time
figure;
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    curNeuron = Simulation.Neuron{iNeuron};
    times = curNeuron.NormalizedData(:,1);
    stimValues = curNeuron.NormalizedData(:,2);
    stimsAfterLinearFilter = curNeuron.NormalizedData(:,3);
    stimsAfterGenerator = curNeuron.NormalizedData(:,4);
    PSTH = curNeuron.PSTH;
    
    %we display only partial data
    dataPoints = 1:(ceil(Simulation.TICKS_IN_SECOND/Simulation.STIMULUS_EACH_TICKS))*SECONDS_OF_RATE_TO_DISPLAY;
    dataPoints = dataPoints';
    times = times(dataPoints);
    stimValues = stimValues(dataPoints);
    stimsAfterLinearFilter = stimsAfterLinearFilter(dataPoints);
    stimsAfterGenerator = stimsAfterGenerator(dataPoints);
       
    %normalize time so we'd start for 0
    onset = times(1);
    times = times-onset;
    PSTH(:,1) = PSTH(:,1)-onset;
    PSTH(PSTH(:,1)<0,:)=[];
    
    %dataPoints(length(dataPoints)+1,1)=length(dataPoints)+1;
    %PSTH = PSTH(dataPoints, :);
    PSTH = PSTH(1:length(dataPoints)-1, :);  
    
    PSTH = [0 0; PSTH];
    
    comp = [times PSTH]; %TODO: remove when stable
    
    timesOfPSTH = PSTH(:,1);
    PSTH = PSTH(:,2); %get values
    NORMALIZE = 1;
    NORMALIZE_BY_MAX = 1;
    if (NORMALIZE) 
        m1 = mean(PSTH);
        
        PSTH = PSTH-m1;
        
        %normalize
        if (NORMALIZE_BY_MAX) 
            PSTH = (PSTH)/max(abs(PSTH));
        end
    end
    
         
    hold on;
           
    plot(times,stimsAfterLinearFilter,'k');
    %h = plot(times,stimValues, 'b');
    h = plot(timesOfPSTH,PSTH, 'b');
    h.Color(4) = 0.3;  % 70% transparent
    h = plot(times,stimsAfterGenerator,'g');
    h.Color(4) = 0.5;  % 50% transparent
        
    title(sprintf('Neuron #%d', iNeuron));
    legend('After STA Linear Filter', ...
    sprintf('PSTH (Bin Size = %.2f ms)', ...
        Simulation.PSTH_BIN_SIZE/Simulation.TICKS_IN_SECOND*1000), ...
    'After Generator (Rest)');
    xlim([0,times(end)]);
    
    %{
    set(gca,'XTickLabel',sprintf('%1.0f|',...
        0:SECONDS_IN_WINDOW:times(end)/TICKS_IN_SECOND));
    %}

    xlabel('Time (s)');
    ylabel('Normalized Values');
    
    hold off;
   
end