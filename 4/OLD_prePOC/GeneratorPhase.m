close all;

ConstantsHeader();

%choose Rep or NonRep
MODE = 'Rep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)) || ...
        Simulation.Phase < CONSTANTS.PHASES.LINEARFILTER)
    clearvars -except MODE;
    load(['AfterLinearFilter_' MODE '.mat'])
    ConstantsHeader();
end

switch MODE
    case 'Rep'
        StimTime = Simulation.StimTimeRep;
    case 'NonRep'
        StimTime = Simulation.StimTimeNonRep;
    otherwise
        ME = MException('STA:noSuchMODE', ...
            'no such MODE is found!');
        throw(ME)
end

SAVE_MAT_FILE = 1;

NEURONS = length(Simulation.Neuron);
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;
SECONDS_OF_RATE_TO_DISPLAY = Simulation.SECONDS_OF_RATE_TO_DISPLAY;

%store for later usage
Simulation.Phase = CONSTANTS.PHASES.LINEARFILTER;


for iNeuron=1:NEURONS
    curNeuron = Simulation.Neuron{iNeuron};
    
    %% create the fit
    data = curNeuron.CurveFitData;
    
    data = sortrows(data, 1);
    XVals = data(:,1); %after linear filter
    YVals = data(:,2); %raw stims

    BIN_SIZE=0.1;

    %create nice looking numbers
    firstBin = floor(min(XVals)*10)/10;
    lastBin = ceil(max(XVals)*10)/10;
    
    %put into bins
    bins = firstBin:BIN_SIZE:lastBin;
    [bincounts,binIndex] = histc(XVals,bins);

    %group by mean
    funcXData = bins';
    funcYMeans = accumarray(binIndex, YVals, [length(bins) 1], @mean);
    
    %remove empty/'zero' intesity buckets
    m = [funcXData funcYMeans];
    m = m(m(:,2)~=0,:);
    funcXData = m(:,1);
    funcYMeans = m(:,2);

    %create the 'function' using Interpolant Linear Fit
    [fitresult, gof] = createFit(funcXData,funcYMeans,iNeuron,MODE);
    curNeuron.Generator = fitresult;
    
    Simulation.Neuron{iNeuron} = curNeuron;    
    
    %% apply generator

    %create column for filtered data
    Simulation.Neuron{iNeuron}.Data(:,5) = ...
        NaN(length(curNeuron.Data(:,1)),1);
        
    stimsAfterNonLinearFilter = curNeuron.Data(:,4);
    res = curNeuron.Generator(stimsAfterNonLinearFilter);
    
    curNeuron.Data(:,5) = res;
    
    
    Simulation.Neuron{iNeuron} = curNeuron;
    
end %for iNeuron

%% normalizing data
for iNeuron=1:NEURONS
    curNeuron = Simulation.Neuron{iNeuron};
    times = curNeuron.Data(:,1);
    stimValues = curNeuron.Data(:,2);
    stimsAfterLinearFilter = curNeuron.Data(:,4);
    stimsAfterGenerator = curNeuron.Data(:,5);
        
    sim_start = StimTime(1);
    sim_end = StimTime(end);
    
    %discard 0-padded conv data
    filter = logical(times(:) >= sim_start & times(:) < sim_end ... 
        & ~isnan(stimValues(:)) ...
        & ~isnan(stimsAfterLinearFilter(:)) ...
        & ~isnan(stimsAfterGenerator(:)));
    windowData = curNeuron.Data(filter, :);
    
    times = windowData(:,1);
    stimValues = windowData(:,2);
    stimsAfterLinearFilter = windowData(:,4);
    stimsAfterGenerator = windowData(:,5);
    
    NORMALIZE = 1;
    NORMALIZE_BY_MAX = 1;
    if (NORMALIZE) 
        m1 = mean(stimValues);
        m2 = mean(stimsAfterLinearFilter);
        m3 = mean(stimsAfterGenerator);

        stimValues = stimValues-m1;
        stimsAfterLinearFilter = stimsAfterLinearFilter-m2;
        stimsAfterGenerator = stimsAfterGenerator-m3;
        
        %normalize
        if (NORMALIZE_BY_MAX) 
            stimValues = (stimValues)/max(abs(stimValues));
            stimsAfterLinearFilter = (stimsAfterLinearFilter)/max(abs(stimsAfterLinearFilter));
            stimsAfterGenerator = (stimsAfterGenerator)/max(abs(stimsAfterGenerator));
        end
    end
    
    curNeuron.NormalizedData = NaN(length(stimValues), 4);
    curNeuron.NormalizedData(:,1) = times;
    curNeuron.NormalizedData(:,2) = stimValues;
    curNeuron.NormalizedData(:,3) = stimsAfterLinearFilter;
    curNeuron.NormalizedData(:,4) = stimsAfterGenerator;
    
    Simulation.Neuron{iNeuron} = curNeuron;
end

%% plotting stims and after filters stims against time
figure;
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    curNeuron = Simulation.Neuron{iNeuron};
    times = curNeuron.NormalizedData(:,1);
    stimValues = curNeuron.NormalizedData(:,2);
    stimsAfterLinearFilter = curNeuron.NormalizedData(:,3);
    stimsAfterGenerator = curNeuron.NormalizedData(:,4);
         
    %we display only partial data
    dataPoints = 1:(ceil(Simulation.TICKS_IN_SECOND/Simulation.STIMULUS_EACH_TICKS))*SECONDS_OF_RATE_TO_DISPLAY;
    times = times(dataPoints);
    stimValues = stimValues(dataPoints);
    stimsAfterLinearFilter = stimsAfterLinearFilter(dataPoints);
    stimsAfterGenerator = stimsAfterGenerator(dataPoints);
    
    %normalize time so we'd start for 0
    times = times-times(1);
         
    hold on;
           
    plot(times,stimsAfterLinearFilter,'k');
    h = plot(times,stimValues, 'b');
    h.Color(4) = 0.3;  % 70% transparent
    h = plot(times,stimsAfterGenerator,'g');
    h.Color(4) = 0.5;  % 50% transparent
        
    title(sprintf('Neuron #%d', iNeuron));
    legend('After STA Linear Filter', 'Raw', 'After Generator');
    xlim([0,times(end)]);
    
    %{
    set(gca,'XTickLabel',sprintf('%1.0f|',...
        0:SECONDS_IN_WINDOW:times(end)/TICKS_IN_SECOND));
    %}

    xlabel('Time (s)');
    ylabel('Normalized Stimuli Values');
    
    hold off;
   
end

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['AfterGenerator_' MODE '.mat'], 'Simulation');
end
        
beep('on');
