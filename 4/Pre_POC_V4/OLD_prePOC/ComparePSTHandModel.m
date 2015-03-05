close all;

ConstantsHeader();

%on question 3, we where asked to compare PSTH to model in Rep only
MODE = 'Rep';
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
BIN_SIZES = Simulation.PSTH_BIN_SIZES;
BIN_SIZE_TO_PLOT = BIN_SIZES(4);

%% plotting stims and after filters stims against time
figure;
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    curNeuron = Simulation.Neuron{iNeuron};
    times = curNeuron.NormalizedData(:,1);
    stimValues = curNeuron.NormalizedData(:,2);
    stimsAfterLinearFilter = curNeuron.NormalizedData(:,3);
    stimsAfterGenerator = curNeuron.NormalizedData(:,4);
    
    %TODO: plot only the binsize to display and for other only store data
    PSTH = curNeuron.PSTH{4};
    
    NORMALIZE = 1;
    NORMALIZE_BY_MAX = 1;
    if (NORMALIZE) 
        m1 = mean(PSTH(:,2));
        
        PSTH(:,2) = PSTH(:,2)-m1;
        
        %normalize
        if (NORMALIZE_BY_MAX) 
            PSTH(:,2) = (PSTH(:,2))/max(abs(PSTH(:,2)));
        end
    end
    
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
    PSTH = PSTH(:,2);
        
         
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
        BIN_SIZE_TO_PLOT/Simulation.TICKS_IN_SECOND*1000), ...
    'After Generator (Rest)');
    xlim([0,times(end)]);
    
    %{
    set(gca,'XTickLabel',sprintf('%1.0f|',...
        0:SECONDS_IN_WINDOW:times(end)/TICKS_IN_SECOND));
    %}

    xlabel('Time (s)');
    ylabel('Normalized Values');
    
    hold off;
    
    %see http://www.mathworks.com/matlabcentral/answers/104189-calculate-sum-of-square-error
    SSEk = norm(PSTH-stimsAfterLinearFilter,2)^2; %lower is less err
    SSEg = norm(PSTH-stimsAfterGenerator,2)^2; %lower is less err

    %the similarity between two signals, we only need zero lag
    %CC = xcorr(PSTH,stimsAfterGenerator,0,'coeff'); % 1 if are equal
    
    text(0, 0.9, ...
        sprintf('\bSSE after STA kernel: %.2f;\nSSE after Generator: %.2f;', ...
            SSEk, SSEg));
   
end