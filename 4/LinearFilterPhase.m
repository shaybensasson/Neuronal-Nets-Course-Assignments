close all;

ConstantsHeader();

%choose Rep or NonRep
MODE = 'Rep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)) || ...
        Simulation.Phase < CONSTANTS.PHASES.STASTC)
    clearvars -except MODE;
    load(['AfterSTA_' MODE '.mat'])
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
ITERATIONS = Simulation.ITERATIONS;
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;

STA_WINDOW_IN_TICKS = Simulation.STA_WINDOW_IN_TICKS;
STIMS_IN_STA_WINDOW = Simulation.STIMS_IN_STA_WINDOW;

%throw 0-padded convolved values, see conv doc for more info;
STIMS_TO_THROW_AFTER_CONV = ceil(STIMS_IN_STA_WINDOW/2);
%duration in seconds to display when ploting multiple rates
SECONDS_OF_RATE_TO_DISPLAY = 10;

%store for later usage
Simulation.Phase = CONSTANTS.PHASES.LINEARFILTER;
Simulation.STIMS_TO_THROW_AFTER_CONV = STIMS_TO_THROW_AFTER_CONV;
Simulation.SECONDS_OF_RATE_TO_DISPLAY = SECONDS_OF_RATE_TO_DISPLAY;
        
%% apply a linear filter of STA
for iNeuron=1:NEURONS
    curNeuron = Simulation.Neuron{iNeuron};
    STA = curNeuron.STA;
    
    %create column for filtered data
    Simulation.Neuron{iNeuron}.Data(:,4) = ...
        NaN(length(Simulation.Neuron{iNeuron}.Data(:,1)),1);
        
    for iIteration=2:ITERATIONS
        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 
 
        times = curNeuron.Data(:,1);
        stimValues = curNeuron.Data(:,2);
        
        onset = StimTime(iStart);
        next_onset = StimTime(iEnd);
        
        %filter window of stimValues
        filter = logical(times(:) >= onset & times(:) < next_onset ...
            & ~isnan(stimValues(:)));
        stimValuesOnWindow = curNeuron.Data(filter, 2);
                
        %{ 
        linear filter the stimVals by the STA filter
        we flip the STA so the convolution will be correct
        
        see http://en.wikipedia.org/wiki/Cross-correlation
        %}
        stimsAfterLinearFilter = conv(stimValuesOnWindow, flipud(STA'), 'valid');

        %added NaNs to 0-padded convolved values, see conv doc for more info;
        %FUTURE: we might have problems here with different sizes of STA
        stimsAfterLinearFilter = padarray(stimsAfterLinearFilter,...
            [STIMS_TO_THROW_AFTER_CONV-1 0], NaN, 'pre');
        stimsAfterLinearFilter = padarray(stimsAfterLinearFilter,...
            [STIMS_TO_THROW_AFTER_CONV 0], NaN, 'post');
        
        curNeuron.Data(filter,4) = stimsAfterLinearFilter;
        
        filter = isnan(curNeuron.Data(:,3)) & curNeuron.Data(:,4)==0;
        curNeuron.Data(filter,4)=NaN; %convert zeros to NaN
        
    end %for ITERATIONS
    
    Simulation.Neuron{iNeuron} = curNeuron;
    
end %for iNeuron

%% plotting raw vs after filter stims and more data analysis
figure;
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    curNeuron = Simulation.Neuron{iNeuron};
    times = curNeuron.Data(:,1);
    stimValues = curNeuron.Data(:,2);
    stimsAfterLinearFilter = curNeuron.Data(:,4);
        
    sim_start = StimTime(1);
    sim_end = StimTime(end);
    
    %discard 0-padded conv data
    filter = logical(times(:) >= sim_start & times(:) < sim_end ... 
        & ~isnan(stimValues(:)) ...
        & ~isnan(stimsAfterLinearFilter(:)));
    windowData = curNeuron.Data(filter, :);
    
    times = windowData(:,1);
    stimValues = windowData(:,2);
    stimsAfterLinearFilter = windowData(:,4);
    
    hold on;
    
    NORMALIZE = 1;
    NORMALIZE_BY_MAX = 1;
    if (NORMALIZE) 
        m1 = mean(stimValues);
        m2 = mean(stimsAfterLinearFilter);

        stimValues = stimValues-m1;
        stimsAfterLinearFilter = stimsAfterLinearFilter-m2;
        
        %normalize
        if (NORMALIZE_BY_MAX) 
            stimValues = (stimValues)/max(abs(stimValues));
            stimsAfterLinearFilter = (stimsAfterLinearFilter)/max(abs(stimsAfterLinearFilter));
        end
    end
    
    curNeuron.NormalizedData = NaN(length(stimValues), 3);
    curNeuron.NormalizedData(:,1) = times;
    curNeuron.NormalizedData(:,2) = stimValues;
    curNeuron.NormalizedData(:,3) = stimsAfterLinearFilter;
    
    
    title(sprintf('Neuron #%d', iNeuron));
    plot(stimsAfterLinearFilter,stimValues,'.');
   
    
    curNeuron.CurveFitData = [stimsAfterLinearFilter stimValues];
    
    Simulation.Neuron{iNeuron} = curNeuron;
    
end

%% plotting stims and after filter stims against time
figure;
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    %figure;
    
    curNeuron = Simulation.Neuron{iNeuron};
    times = curNeuron.NormalizedData(:,1);
    stimValues = curNeuron.NormalizedData(:,2);
    stimsAfterLinearFilter = curNeuron.NormalizedData(:,3);
    
    %we display only partial data
    dataPoints = 1:(ceil(Simulation.TICKS_IN_SECOND/Simulation.STIMULUS_EACH_TICKS))*SECONDS_OF_RATE_TO_DISPLAY;
    times = times(dataPoints);
    stimValues = stimValues(dataPoints);
    stimsAfterLinearFilter = stimsAfterLinearFilter(dataPoints);
    
    %normalize time so we'd start for 0
    times = times-times(1);
         
    hold on;
           
    plot(times,stimsAfterLinearFilter,'k');
    h = plot(times,stimValues, 'b');
    h.Color(4) = 0.3;  % 50% transparent
        
    title(sprintf('Neuron #%d', iNeuron));
    legend('After STA Linear Filter', 'Raw');
    
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
    save(['AfterLinearFilter_' MODE '.mat'], 'Simulation');
end
        
beep('on');

