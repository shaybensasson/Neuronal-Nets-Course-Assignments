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
ITERATIONS = Simulation.ITERATIONS;
PSTH_BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq
             Simulation.STIMULUS_EACH_TICKS*2; ... ~66 msec          
             Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec          
             Simulation.STIMULUS_EACH_TICKS*4; ... ~120 msec
             Simulation.STIMULUS_EACH_TICKS*5]; ... ~150 sec
             
STIMULI_PER_WINDOW = Simulation.STIMULI_PER_WINDOW;
%store for later usage
Simulation.Phase = CONSTANTS.PHASES.LINEARFILTER;
Simulation.PSTH_BIN_SIZES = PSTH_BIN_SIZES;

for iBinSize=1:numel(PSTH_BIN_SIZES)
    curBinSize = PSTH_BIN_SIZES(iBinSize);
    
    figure( 'Name', sprintf('Generator (%s), Bin size: %.2f', ...
        MODE, curBinSize));

    for iNeuron=1:NEURONS
        subplot(2,2,iNeuron);
        curNeuron = Simulation.Neuron{iNeuron};
        fprintf('[N:#%i] ...\n', iNeuron);

        idxAcc = 0;
        accBinned = NaN(STIMULI_PER_WINDOW*ITERATIONS,2);
        for iIteration=2:ITERATIONS
            fprintf('[N:#%i] processing iteration #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

            %Normalize stim time to start from 1
            iEnd = iIteration;
            iStart = iEnd-1; 

            times = Simulation.Neuron{iNeuron}.Data(:,1);

            onset = StimTime(iStart);
            next_onset = StimTime(iEnd);
            filter = logical(times(:) >= onset & times(:) < next_onset);
            windowData = curNeuron.Data(filter, :);

            times = windowData(:, 1);
            times = times - onset+1;
            stimValues = windowData(:,2);
            aps = windowData(:,3);
            stimsAfterLinearFilter = windowData(:,4);

            %bins = times(start):BIN_SIZE:times(end);
            bins = 0:curBinSize:TICKS_IN_WINDOW;
            [~,binIndex] = histc(times(:,1),bins);
            %insert any unbinned data to last bin
            binIndex(binIndex==0)=max(binIndex);

            %group by mean and sum
            sumIgnoreNaNs = @(vector) sum(vector(~isnan(vector(:))));
            meanIgnoreNaNs = @(vector) mean(vector(~isnan(vector(:))));

            %grp_stimValues = accumarray(binIndex, stimValues, [length(bins) 1], meanIgnoreNaNs);
            grp_aps = accumarray(binIndex, aps, [length(bins) 1], sumIgnoreNaNs);
            grp_stimsAfterLinearFilter = accumarray(binIndex, stimsAfterLinearFilter, [length(bins) 1], meanIgnoreNaNs);

            %binned = [ceil(bins') grp_stimValues grp_aps grp_stimsAfterLinearFilter];
            binned = [grp_stimsAfterLinearFilter grp_aps];
            binned(binned(:,1) == 0,:) = [];
            accBinned(idxAcc+1:idxAcc+length(binned),:) = binned;
            idxAcc = idxAcc+length(binned);
        end

        %remove unused preallocated rows
        accBinned(accBinned(:,1) == 0 & accBinned(:,2) == 0, :) = [];

        %% normalize
        stimsAfterLinearFilter = accBinned(:,1);
        psth = accBinned(:,2);

        NORMALIZE = 1;
        NORMALIZE_BY_MAX = 1;
        stimsAfterLinearFilter = normalize(stimsAfterLinearFilter, NORMALIZE, NORMALIZE_BY_MAX);
        psth = normalize(psth, NORMALIZE, NORMALIZE_BY_MAX);

        %% create the fit
        data = [stimsAfterLinearFilter psth];
        data(isnan(data(:,1)) | isnan(data(:,2)), :) = [];

        data = sortrows(data, 1);
        XVals = data(:,1); %after linear filter
        YVals = data(:,2); %psth

        FIT_BIN_SIZE=0.1;

        %create nice looking numbers
        firstBin = floor(min(XVals)*10)/10;
        lastBin = ceil(max(XVals)*10)/10;

        %put into bins
        bins = firstBin:FIT_BIN_SIZE:lastBin;
        [bincounts,binIndex] = histc(XVals,bins);

        %group by mean
        funcXData = bins';
        funcYMeans = accumarray(binIndex, YVals, [length(bins) 1], @mean);

        %remove empty/'zero' intesity buckets
        m = [funcXData funcYMeans];
        %m = m(m(:,2)~=0,:);
        funcXData = m(:,1);
        funcYMeans = m(:,2);

        %create the 'function' using Interpolant Linear Fit
        [fitresult, gof] = createFit(funcXData,funcYMeans,iNeuron,MODE);
        curNeuron.PSTH{iBinSize}.BIN_SIZE = curBinSize;
        curNeuron.PSTH{iBinSize}.Generator = fitresult;
                
        Simulation.Neuron{iNeuron} = curNeuron;    
    end %for iNeuron
    
    %% apply the generator
    for iNeuron=1:NEURONS
        curNeuron = Simulation.Neuron{iNeuron};
        fprintf('[N:#%i] ...\n', iNeuron);
        
        stimsAfterNonLinearFilter = curNeuron.Data(:,4);
        NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
        stimsAfterNonLinearFilter = normalize(stimsAfterNonLinearFilter, NORMALIZE, NORMALIZE_BY_MAX);
        res = curNeuron.PSTH{iBinSize}.Generator(stimsAfterNonLinearFilter);

        curNeuron.PSTH{iBinSize}.Result = res; 
        
        Simulation.Neuron{iNeuron} = curNeuron;
    end
end %for iBinSize

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['AfterGenerator_' MODE '.mat'], 'Simulation');
end
        
beep('on');
