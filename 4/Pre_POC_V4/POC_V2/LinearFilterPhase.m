close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

load(['AfterSTA_' MODE '.mat'])
    
switch MODE
    case 'Rep'
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
    fprintf('[N:#%i] ...\n', iNeuron);
    
    %create Data that contains normalized data and a column for the filtered data
    curNeuron = Simulation.Neuron{iNeuron};
    curNeuron.Data = Simulation.Neuron{iNeuron}.RawData;
    curNeuron.Data(:, 4) = NaN(length(curNeuron.Data), 1);

    STA = curNeuron.STA;

    times = curNeuron.RawData(:,1);
    stimValues = curNeuron.RawData(:,2);
    
    NORMALIZE = 1; NORMALIZE_BY_MAX=1;
    stimValues = normalize(stimValues, NORMALIZE, NORMALIZE_BY_MAX);
    %update normalized values
    curNeuron.Data(:,2) = stimValues;
        
    for iIteration=2:ITERATIONS
        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 
         
        onset = StimTime(iStart);
        next_onset = StimTime(iEnd);
        
        %filter window of stimValues
        filter = logical(times(:) >= onset & times(:) < next_onset ...
            & ~isnan(stimValues(:)));
        stimValuesOnWindow = curNeuron.RawData(filter, 2);
                
        %{ 
        linear filter the stimVals by the STA filter
        we flip the STA so the convolution will be act as moving window
        
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
    end %for ITERATIONS
    
    %NOTE: somehow the filter inserts zeros on data that 
    %  is not included in the query, so we replace them with NaNs
    %TODO: it might not appear anymore, can be removed?
    %curNeuron.Data(curNeuron.Data(:,4) == 0,4) = NaN;

    %remove conv padded rows - ap column will be empty, and so is the conv
    filter = (isnan(curNeuron.Data(:,3)) & isnan(curNeuron.Data(:,4)));
    curNeuron.Data(filter,:)=[]; %remove padded
        
    Simulation.Neuron{iNeuron} = curNeuron;
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
        grp_psth = accumarray(binIndex, aps, [length(bins) 1], @nansum, NaN);

        %group by times and mean stim vals
        grp_afterLinearFilter = accumarray(binIndex, stimsAfterLinearFilter, [length(bins) 1], @nanmean, NaN);

        rateData = [bins' grp_psth grp_afterLinearFilter ];

        %remove NaNs
        rateData(isnan(rateData(:, 2)) | isnan(rateData(:, 3)), :) = [];

        %throw bins before first ap
        firstApTime = curNeuron.Data(~isnan(aps) & aps>0, 1);
        firstApTime = firstApTime(1);

        %throw bins before first ap
        %FUTURE: why do we check for psth == 0?
        rateData(rateData(:,1)<firstApTime & rateData(:,2) == 0, :) = [];

        %throw last bin, it has garbage unbinned data (it's acc is NaN)
        rateData = rateData(1:end-1,:); 
        psth = rateData(:, 2);
        binnedStimsAfterLinearFilter = rateData(:, 3);

        NORMALIZE = 1;
        NORMALIZE_BY_MAX = 1;

        psth = normalize(psth, NORMALIZE, NORMALIZE_BY_MAX);
        rateData(:, 2) = psth;
        binnedStimsAfterLinearFilter = normalize(binnedStimsAfterLinearFilter, NORMALIZE, NORMALIZE_BY_MAX);
        rateData(:, 3) = binnedStimsAfterLinearFilter;

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
    save(['AfterLinearFilter_' MODE '.mat'], ['Sim_' MODE]);
end
        
beep('on');

