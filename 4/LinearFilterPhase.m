close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';
UsingSTA = 1;

load(['MatFiles\AfterSTA_' MODE '.mat'])
    
switch MODE
    case 'Rep'
        load('MatFiles\AfterSTA_NonRep.mat') %we're using its STA
        Simulation = Sim_Rep;
        clearvars Sim_Rep;
        StimTime = Simulation.StimTimeRep;
    case 'NonRep'
        Simulation = Sim_NonRep;
        clearvars Sim_NonRep;
        StimTime = Simulation.StimTimeNonRep;
    otherwise
        ME = MException('noSuchMODE', ...
            'no such MODE is found!');
        throw(ME)
end

SAVE_MAT_FILE = 1;

NEURONS = length(Simulation.Neuron);
ITERATIONS = Simulation.ITERATIONS;
SECONDS_IN_TRAIL = Simulation.SECONDS_IN_TRAIL;
TICKS_IN_TRAIL = Simulation.TICKS_IN_TRAIL;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;

STA_WINDOW_IN_TICKS = Simulation.STA_WINDOW_IN_TICKS;

%store for later usage
Simulation.Phase = CONSTANTS.PHASES.LINEARFILTER;
Simulation.RATE_BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq
             Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec          
             Simulation.STIMULUS_EACH_TICKS*5; ... ~150 msec
         ];
     
%Determines whether to use STA or STC filters
Simulation.UsingSTA = UsingSTA;
        
%% apply a linear filter of STA
for iNeuron=1:NEURONS
    fprintf('[N:#%i] ...\n', iNeuron);
    
    curNeuron = Simulation.Neuron{iNeuron};
        
    %% get the filter
    if (Simulation.UsingSTA)
        switch MODE
            case 'Rep'
                %IMPORTANT: apply the non-rep STA
                Filter = Sim_NonRep.Neuron{iNeuron}.STA;
            case 'NonRep'
                Filter = Simulation.Neuron{iNeuron}.STA;
        end
    else
        switch MODE
            case 'Rep'
                %IMPORTANT: apply the non-rep STC
                Filter = Sim_NonRep.Neuron{iNeuron}.STCFilter;
            case 'NonRep'
                Filter = Simulation.Neuron{iNeuron}.STCFilter;
        end
    end
    
    %% apply the filter
    numOfTicksInTrail = length(Simulation.Neuron{iNeuron}.Iteration{1});    
    neuronData = NaN(numOfTicksInTrail*ITERATIONS, 4);
    lastIterationIndex = 0;
    
    for iIteration=1:ITERATIONS
        %merge filtered column into data
        data = Simulation.Neuron{iNeuron}.Iteration{iIteration};
        stimValues = data(:,2);
    
        %NOTE: The Generator is indifferent to the normalize of stimValues
        %NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
        %stimValues = normalize(stimValues, NORMALIZE, NORMALIZE_BY_MAX);
        
        % convolve flipped filter (because conv flips it) with stimValues
        stimsAfterLinearFilter = conv(stimValues, flipud(Filter'), 'valid');
        %stimsAfterLinearFilter = conv(stimValues, fliplr(Filter), 'full'); %for full conv
        
        %{
        % run a Filter as a linear window on stimValues
        %stimsAfterLinearFilter2 = funLinearWindow(stimValues, Filter');
        %}

        %{ 
                *********** OLD ************** 
        keeping for lookup and understanding of the whole processs

        % convolve filter with stimValues
        stimsAfterLinearFilter = conv(stimValues, Filter');

        %throw partial convolved vals
        STIMS_TO_THROW_AFTER_CONV = STA_WINDOW_IN_TICKS/2; %5000
        stimsAfterLinearFilter(1:STIMS_TO_THROW_AFTER_CONV-1) = [];
        stimsAfterLinearFilter(end-STIMS_TO_THROW_AFTER_CONV+1:end) = [];
        
        %throw partial convolved vals, throw entire STA window from start
        %stimsAfterLinearFilter(1:STA_WINDOW_IN_TICKS-1) = [];
        %}
        
        %convedDiff = length(stimsAfterLinearFilter)-length(stimValues); %for full conv
        convedDiff = length(stimValues) - length(stimsAfterLinearFilter);
        
        data = data(convedDiff+1:end, :);
        
        %NOTE: normalize stim values fater conv
        NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
        stimsAfterLinearFilter = normalize(stimsAfterLinearFilter, NORMALIZE, NORMALIZE_BY_MAX);

        data(:,4) = stimsAfterLinearFilter;
        %data(end-convedDiff: end, :) = []; %for full conv
        
        %NOTE: free space
        Simulation.Neuron{iNeuron}.Iteration{iIteration} = [];
                
        %store in neuronData output
        %nextIterationIndex = lastIterationIndex + numOfTicksInTrail + 1;
        nextIterationIndex = lastIterationIndex + length(data) + 1;
        neuronData(lastIterationIndex+1:nextIterationIndex-1,:) = data;
        lastIterationIndex = nextIterationIndex-1;
    end %for ITERATIONS
    
    %get only used data vs allocated
    neuronData = neuronData(1:lastIterationIndex, :);

    Simulation.Neuron{iNeuron}.Data = neuronData;
end %for iNeuron

%% binnify into psth bins
fprintf('\nBinnify into psth bins ... \n')
for iBinSize=1:numel(Simulation.RATE_BIN_SIZES)
    curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);

    for iNeuron = 1:NEURONS
        fprintf('[Bsz,N:#%d,#%d] ...\n', iBinSize, iNeuron);
        
        curNeuron = Simulation.Neuron{iNeuron};
        
        times = curNeuron.Data(:,1);
        stimValues = curNeuron.Data(:,2);
        aps = curNeuron.Data(:,3);
        stimsAfterLinearFilter = curNeuron.Data(:,4);

        %group raw data times by bins
        bins = times(1):curBinSize:times(end);
        [~,binIndex] = histc(times,bins);

        %insert any unbinned data to last bin
        binIndex(binIndex==0)=max(binIndex);

        %group by times and sum aps
        %   we put 0 in bins without aps
        grp_psth = accumarray(binIndex, aps, [length(bins) 1], @nansum, 0);

        %group by times and mean stim vals
        %   we put NaN in empty bins (data that was wiped out of iteration)
        grp_afterLinearFilter = accumarray(binIndex, stimsAfterLinearFilter, [length(bins) 1], @mean, NaN);

        rateData = [bins' grp_psth grp_afterLinearFilter];

        %remove NaNs
        rateData(isnan(rateData(:, 3)), :) = [];

        %throw last bin, it has garbage unbinned data (it's acc is NaN)
        rateData = rateData(1:end-1,:); 
        psth = rateData(:, 2);
        binnedStimsAfterLinearFilter = rateData(:, 3);

        %Normalize PSTH and stims after conv
        NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
        psth = normalize(psth, NORMALIZE, NORMALIZE_BY_MAX);
                
        NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
        binnedStimsAfterLinearFilter = normalize(binnedStimsAfterLinearFilter, NORMALIZE, NORMALIZE_BY_MAX);
        
        curNeuron.Rate{iBinSize}.BinSize = curBinSize;
        curNeuron.Rate{iBinSize}.Data = [rateData(:,1) psth binnedStimsAfterLinearFilter];
        
        Simulation.Neuron{iNeuron} = curNeuron;
    end
end

switch MODE
    case 'Rep'
        Sim_Rep = Simulation;
    case 'NonRep'
        Sim_NonRep = Simulation;
    otherwise
        ME = MException('noSuchMODE', ...
            'no such MODE is found!');
        throw(ME)
end

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['MatFiles\AfterLinearFilter_' MODE '.mat'], ['Sim_' MODE], '-v7.3');
end


load gong 
sound(y,Fs)
