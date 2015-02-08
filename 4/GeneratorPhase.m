close all;

ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

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

SAVE_MAT_FILE = 0;

NEURONS = length(Simulation.Neuron);
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;
SECONDS_OF_RATE_TO_DISPLAY = Simulation.SECONDS_OF_RATE_TO_DISPLAY;
ITERATIONS = Simulation.ITERATIONS;
PSTH_BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq (1)
             Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec (2)          
             Simulation.STIMULUS_EACH_TICKS*5; ... ~150 msec (3)
             Simulation.STIMULUS_EACH_TICKS*10; ... ~330 msec (4)
             Simulation.STIMULUS_EACH_TICKS*15; ... ~500 msec (5)
             Simulation.STIMULUS_EACH_TICKS*30]; %1000 msec = 1 sec (6)

             
STIMULI_PER_WINDOW = Simulation.STIMULI_PER_WINDOW;
%store for later usage
Simulation.Phase = CONSTANTS.PHASES.GENERATOR;
Simulation.PSTH_BIN_SIZES = PSTH_BIN_SIZES;

for iBinSize=4:4
%TODO: for iBinSize=1:numel(PSTH_BIN_SIZES)
    curBinSize = PSTH_BIN_SIZES(iBinSize);
    
    figure( 'Name', sprintf('Generator (%s), Bin size: %.2f', ...
        MODE, curBinSize));
    
    for iNeuron=2:2
    %for iNeuron=1:NEURONS
        %TODO: subplot(2,2,iNeuron);
        curNeuron = Simulation.Neuron{iNeuron};
        fprintf('[Bsz,N:#%d,#%d] ...\n', iBinSize, iNeuron);

        %get binned psth and stims after K
        psth = binned(:,2);
        stimsAfterLinearFilter = binned(:,3);
        
        %% create the fit
        curveFitData = [stimsAfterLinearFilter psth];
        curveFitData(isnan(curveFitData(:,1)) | isnan(curveFitData(:,2)), :) = [];

        curveFitData = sortrows(curveFitData, 1);
        XVals = curveFitData(:,1); %after linear filter
        YVals = curveFitData(:,2); %psth

        FIT_BIN_SIZE=0.1;

        %create nice looking numbers
        firstBin = floor(min(XVals)*10)/10;
        lastBin = ceil(max(XVals)*10)/10;

        %put into bins
        bins = firstBin:FIT_BIN_SIZE:lastBin;
        [bincounts,binIndex] = histc(XVals,bins);

        %group by mean
        funcXData = bins';
        funcYMeans = accumarray(binIndex, YVals, [length(bins) 1], @mean, NaN);

        %remove empty/'zero' intesity buckets
        m = [funcXData funcYMeans];
        m = m(~isnan(m(:,2)),:);
        funcXData = m(:,1);
        funcYMeans = m(:,2);

        %create the 'function' using Interpolant Linear Fit
        [fitresult, gof] = createFit(funcXData,funcYMeans,iNeuron,MODE);
        curNeuron.PSTH{iBinSize}.BIN_SIZE = curBinSize;
        curNeuron.PSTH{iBinSize}.Generator = fitresult;
                
        %% apply the generator
        stimsAfterGenerator = normalize( ...
            fitresult(stimsAfterLinearFilter), ...
            NORMALIZE, NORMALIZE_BY_MAX);
        
        %TODO: curNeuron.PSTH{iBinSize}.RVsRest = accBinned;
        
        %% create statistics
                
        %see http://www.mathworks.com/matlabcentral/answers/104189-calculate-sum-of-square-error
        SSEk = norm(psth-stimsAfterLinearFilter,2)^2; %lower is less err
        SSEg = norm(psth-stimsAfterGenerator,2)^2; %lower is less err
                
        curNeuron.PSTH{iBinSize}.SSEk = SSEk;
        curNeuron.PSTH{iBinSize}.SSEg = SSEg;
        
        %TODO: Simulation.Neuron{iNeuron} = curNeuron;    
    end %for iNeuron
end %for iBinSize

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['AfterGenerator_' MODE '.mat'], 'Simulation');
end
        
beep('on');
