clearvars -except STA STCFilter
close all;

MODE = 'NonRep'; %TODO: if changing this, need to review code

%% for quick and dirty
ITERATIONS=1;

%% Pre processing
if (~exist('TTNonRep','var'))
    load('FixedData.mat')
end

g_iBinSize = 1;
%g_BinSize = 500;
g_BinSize = 333;


%secs in stimulus non-rep trail
SECONDS_IN_TRAIL = 100;
TICKS_IN_SECOND = 10000;
TICKS_IN_TRAIL = TICKS_IN_SECOND * SECONDS_IN_TRAIL;
STIMULI_PER_SECOND = 1/30;
%how many ticks there are between stimuli.
STIMULUS_EACH_TICKS = TICKS_IN_SECOND*STIMULI_PER_SECOND; %333.333
STIMULI_PER_TRAIL = SECONDS_IN_TRAIL / STIMULI_PER_SECOND; %3000 stims for 100 secs

%throw 30 secs on the begining of each iteration, malformed data
TICKS_TO_THROW = 30 * TICKS_IN_SECOND;

iNeuron = 2;

fprintf('[N:#%i] ...\n', iNeuron);

lastIterationSimuliIndex = 0;    

iIteration = 1;
%fprintf('[N:#%i] iter #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

iStart = iIteration;
iEnd = iIteration+1;

%actual ticks in trail limited by TICKS_IN_TRAIL
dataTicks = TICKS_IN_TRAIL;

iteration.time= (StimTimeNonRep(iStart):StimTimeNonRep(iEnd)-1)'; %first 100 secs
dataTicks = min(dataTicks,length(iteration.time)); 

% get stim values of the current non-rep trail        
nextIterationSimuliIndex = lastIterationSimuliIndex + STIMULI_PER_TRAIL + 1;
stimValues = StimulusNonRep(lastIterationSimuliIndex+1:nextIterationSimuliIndex-1,:);
lastIterationSimuliIndex = nextIterationSimuliIndex-1;

%Smoothen stim values on dataTicks
stimValues = repmat(stimValues,1,floor(STIMULUS_EACH_TICKS));
stimValues = stimValues';
stimValues = stimValues(:);

%because of flooring we might have less data
dataTicks = min(dataTicks,length(stimValues)); 

%trim data collected so far
iteration.time = iteration.time(1:dataTicks);
iteration.stimuli=stimValues(1:dataTicks);

% get aps of the current non-rep trail        
timeOfAPs=TTNonRep(iNeuron).sp; %time of Aps

onset = StimTimeNonRep(iStart);
next_onset = StimTimeNonRep(iEnd);
%{
        Change from time scale into discrete indexes and 
        filter by AP occurences

        logical() is required to transform binary vector to a valid filter.
        %}
indexes = 1:length(timeOfAPs);
filter = logical(timeOfAPs(:) >= onset & timeOfAPs(:) < next_onset);
indexesOfAPs = indexes(filter);

iteration.APs = NaN(dataTicks,1);
%convert times to indexes
indexesOfAPsInTrail = timeOfAPs(indexesOfAPs)-onset + 1; %the index is 1 based
indexesOfAPsInTrail(indexesOfAPsInTrail>dataTicks)=[]; %remove aps out of ticks range
iteration.APs(indexesOfAPsInTrail)=1;

%create fixed data set
iteration.all = NaN(dataTicks, 3);
iteration.all(:,1) = iteration.time;
iteration.all(:,2) = iteration.stimuli;
iteration.all(:,3) = iteration.APs;

%throw initial ticks
iteration.all(1:TICKS_TO_THROW,:) = [];

Simulation.Neuron{iNeuron}.Iteration{iStart} = iteration.all;


%% STA
STA_WINDOW_IN_MS = 1000;
STA_WINDOW_IN_SEC = STA_WINDOW_IN_MS/1000;
STA_WINDOW_IN_TICKS = STA_WINDOW_IN_SEC*TICKS_IN_SECOND;

if (~exist('STA','var'))
    load('MatFiles\AfterSTA_NonRep.mat'); %13 secs time
    
    STA = Sim_NonRep.Neuron{iNeuron}.STA;
    STCFilter = Sim_NonRep.Neuron{iNeuron}.STCFilter;
    clearvars Sim_NonRep;
end
Simulation.Neuron{iNeuron}.STA = STA;
Simulation.Neuron{iNeuron}.STCFilter = STCFilter;

%% Linear Filter
UsingSTA = 1;

%{
Simulation.RATE_BIN_SIZES = [STIMULUS_EACH_TICKS; ... %sampling freq
             STIMULUS_EACH_TICKS*3; ... 100 msec          
             STIMULUS_EACH_TICKS*5; ... ~150 msec
         ];
%}

%Determines whether to use STA or STC filters
Simulation.UsingSTA = UsingSTA;

% apply a linear filter of STA
iNeuron = 2;

fprintf('[N:#%i] ...\n', iNeuron);

curNeuron = Simulation.Neuron{iNeuron};

% get the filter
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

%NOTE: normalize the filter
NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
Filter = normalize(Filter, NORMALIZE, NORMALIZE_BY_MAX);

% apply the filter
numOfTicksInTrail = length(Simulation.Neuron{iNeuron}.Iteration{1});    
%neuronData = NaN(numOfTicksInTrail*ITERATIONS, 4);
neuronData = NaN(numOfTicksInTrail*1, 4);
lastIterationIndex = 0;

iIteration=1;

%create another column for the filtered data
data = Simulation.Neuron{iNeuron}.Iteration{iIteration};
stimValues = data(:, 2);

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



%convedDiff = length(stimsAfterLinearFilter)-length(stimValues); %for full conv
convedDiff = length(stimValues) - length(stimsAfterLinearFilter);

data = data(convedDiff+1:end, :);

%NOTE: normalize stim values fater conv
NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
stimsAfterLinearFilter = normalize(stimsAfterLinearFilter, NORMALIZE, NORMALIZE_BY_MAX);

data(:,4) = stimsAfterLinearFilter;
%data(end-convedDiff: end, :) = []; %for full conv

%NOTE: free space
%Simulation.Neuron{iNeuron}.Iteration{iIteration} = [];

%store in neuronData output
%nextIterationIndex = lastIterationIndex + numOfTicksInTrail + 1;
nextIterationIndex = lastIterationIndex + length(data) + 1;
neuronData(lastIterationIndex+1:nextIterationIndex-1,:) = data;
lastIterationIndex = nextIterationIndex-1;

%get only used data vs allocated
neuronData = neuronData(1:lastIterationIndex, :);
    
Simulation.Neuron{iNeuron}.Data = neuronData;

% binnify into psth bins
fprintf('\nBinnify into psth bins ... \n')
%iBinSize=1;
iBinSize=g_iBinSize;
    %curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    curBinSize = g_BinSize;

    iNeuron = 2;
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
        %grp_afterLinearFilter = accumarray(binIndex, stimsAfterLinearFilter, [length(bins) 1], @sum, NaN);

        rateData = [bins' grp_psth grp_afterLinearFilter];

        %remove NaNs
        rateData(isnan(rateData(:, 3)), :) = [];

        %throw last bin, it has garbage unbinned data (it's acc is NaN)
        rateData = rateData(1:end-1,:); 
        times = rateData(:, 1);
        psth = rateData(:, 2);
        binnedStimsAfterLinearFilter = rateData(:, 3);

        %Normalize PSTH and stims after conv
        NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
        psth = normalize(psth, NORMALIZE, NORMALIZE_BY_MAX);
               
        NORMALIZE = 1; NORMALIZE_BY_MAX = 1;
        binnedStimsAfterLinearFilter = normalize(binnedStimsAfterLinearFilter, NORMALIZE, NORMALIZE_BY_MAX);
        
        %TODO: we try abs here to make the estimate positive
        %binnedStimsAfterLinearFilter = abs(binnedStimsAfterLinearFilter);
        
        curNeuron.Rate{iBinSize}.BinSize = curBinSize;
        curNeuron.Rate{iBinSize}.Data = [rateData(:,1) psth binnedStimsAfterLinearFilter];
        
        Simulation.Neuron{iNeuron} = curNeuron;

%% Generator
%iBinSize=1;
iBinSize=g_iBinSize;

    %curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    curBinSize = g_BinSize;
    
    title = sprintf('Generator (%s), Bin size: %.2f, Using %s Filter', ...
        MODE, curBinSize, iif(Simulation.UsingSTA,'STA','STC'));
    
    %TODO: hf = figure(iBinSize);
    hf = figure;
    hf.Name = title;
    
    iNeuron=2;
        %TODO: subplot(2,2,iNeuron);
        curNeuron = Simulation.Neuron{iNeuron};
        fprintf('[Bsz,N:#%d,#%d] ...\n', iBinSize, iNeuron);

        %get binned psth and stims after K
        rateData = curNeuron.Rate{iBinSize}.Data;
        psth = rateData(:,2);
        stimsAfterLinearFilter = rateData(:,3);
        
        % create the fit
        curveFitData = [stimsAfterLinearFilter psth];

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
        
        %insert any unbinned data to last bin (we're gonna remove it later)
        %binIndex(binIndex==0)=max(binIndex);

        %group by mean
        funcXData = bins';
        funcYMeans = accumarray(binIndex, YVals, [length(bins) 1], @mean, NaN);

        %remove empty/'zero' intesity buckets
        m = [funcXData funcYMeans bincounts];
        %m = m(2:end-1, :); %throw first and last bins, really few stims there
        m = m(~isnan(m(:,2)),:);
        
               
        %smoothen the edges with values
        idxs = 1:length(m);
        m = [m idxs'];
        filter = m(:,1)>=-0.6 & m(:,1)<=0.6;
        mFiltered = m(filter, :)
        idxFrom = mFiltered(1,4);
        m(1:idxFrom, 2) = m(idxFrom, 2);
        idxTo = mFiltered(end,4);
        m(idxTo:end, 2) = m(idxTo, 2);

        funcXData = m(:,1);
        funcYMeans = m(:,2);
        bincounts = m(:,3);
                
        %create the 'function' using Interpolant Linear Fit
        [fitresult, gof] = createFit(funcXData,funcYMeans,bincounts,iNeuron,MODE);
        curNeuron.Rate{iBinSize}.Generator = fitresult;
                
        
        % apply the generator
        NORMALIZE = 1; NORMALIZE_BY_MAX=1;
        stimsAfterGenerator = normalize( ...
            fitresult(stimsAfterLinearFilter), ...
            NORMALIZE, NORMALIZE_BY_MAX);
        
        rateData(:,4) = stimsAfterGenerator;
        curNeuron.Rate{iBinSize}.Data = rateData;
        
        % create statistics
                
        %see http://www.mathworks.com/matlabcentral/answers/104189-calculate-sum-of-square-error
        SSEk = norm(psth-stimsAfterLinearFilter,2)^2; %lower is less err
        SSEg = norm(psth-stimsAfterGenerator,2)^2; %lower is less err
                
        curNeuron.Rate{iBinSize}.SSEk = SSEk;
        curNeuron.Rate{iBinSize}.SSEg = SSEg;
        
        Simulation.Neuron{iNeuron} = curNeuron;    
    
    CreateTitleForSubplots(['\bf ' title]);
    
    r = 150; %pixels pre inch
    set(hf, 'PaperUnits', 'inches');
    set(hf, 'PaperPosition', [0 0 2880 1620]/r); %x_width=10cm y_width=15cm
    
    saveas(hf, ['QnD_Generator_' MODE '_BinSize_' ...
        sprintf('%d', floor(curBinSize)) ...
        '_UsingSTA_' num2str(Simulation.UsingSTA)], 'png');
    

%% plot after generator
%duration in seconds to display when ploting multiple rates
LENGTH_OF_PLOT_IN_SECS = 10;
LENGTH_OF_PLOT_IN_TICKS = LENGTH_OF_PLOT_IN_SECS*TICKS_IN_SECOND; %~ num of records
START_PLOT_FROM_INDEX = 1;

COLOR_MAP = [43, 87, 154; ... %blue: 1
            32, 162, 58; ... %green: 2
            0, 0, 0; ... % black: 3
            175, 185, 22; ... %yello for stims: 4
            201, 67, 67]; %red for aps: 5

COLOR_MAP = COLOR_MAP./255;



%iBinSize=1;
iBinSize=g_iBinSize;
    %curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    curBinSize = g_BinSize;
    
    title = sprintf('R vs R_{Est} (%s), Bin size: %.2f, Using %s Filter', ...
        MODE, curBinSize, iif(Simulation.UsingSTA,'STA','STC'));
    
    %TODO: hf = figure(iBinSize);
    hf = figure;
    hf.Name = title;
    
    iNeuron=2;
    %for iNeuron=1:NEURONS
        %TODO: hs = subplot(2,2,iNeuron);
        
        %TODO: PositionOfSubplot = hs.Position;
        PositionOfSubplot = get(gca, 'Position');
        
        curNeuron = Simulation.Neuron{iNeuron};

        data = curNeuron.Data;
        rateData = curNeuron.Rate{iBinSize}.Data;
        
        %% plot
        indexesToPlot = START_PLOT_FROM_INDEX:START_PLOT_FROM_INDEX+LENGTH_OF_PLOT_IN_TICKS-1;
        dataToPlot = data(indexesToPlot, :);
        times = dataToPlot(:,1);
        stimValues = dataToPlot(:,2);
        aps = dataToPlot(:,3);

        xRateFilter = rateData(:,1)>=times(START_PLOT_FROM_INDEX) & ...
            rateData(:,1)<=times(START_PLOT_FROM_INDEX)+LENGTH_OF_PLOT_IN_TICKS;
        binnedToPlot = rateData(xRateFilter, :);
        
        psth = binnedToPlot(:,2);
        stimsAfterLinearFilter = binnedToPlot(:,3);
        stimsAfterGenerator = binnedToPlot(:,4);

        
        %bind aps to top
        minValue = min(min(binnedToPlot(:,2:4)));
        maxValue = max(max(binnedToPlot(:,2:4)));
        
        aps = aps*maxValue;

        hold on;
        [AX,H1,H2] = plotyy(times, stimValues, times, aps, 'plot', 'plot');
        %hold(AX(1));
        set(AX,'NextPlot','add')
        H1.LineStyle = 'none'; H1.Marker = '.'; 
        H1.Color = COLOR_MAP(4, :);  
        
        AX(1).YLim = [min(stimValues) max(stimValues)];
        AX(1).YColor = COLOR_MAP(4, :);
        %AX(1).YTick = linspace(min(stimValues), max(stimValues), 5);

        H2.LineStyle = 'none'; H2.Marker = 'o';
        H2.Color = COLOR_MAP(5, :);
        
        AX(2).YLim = [minValue maxValue];
        AX(2).YColor = COLOR_MAP(1, :);
        %AX(2).YTick = linspace(minValue, maxValue, 5);
        

        %filter relevant data to plot
        AX(1).XLim = [times(1) times(end)];
        AX(2).XLim = [times(1) times(end)];
        
        %{
        ax = gca;
        ticks = ax.XTick;
        
        set(gca,'XTickLabel',sprintf('%.1f\n',...
            ticks ...
                /Simulation.TICKS_IN_SECOND));
        %}
        
        %psth
        h = plot(AX(2), binnedToPlot(:,1), psth);
        h.Color = COLOR_MAP(1, :);
        h.Color(4) = 0.60; % 30% transparent 
        h.LineWidth = 1.5;
        %K filter
        h = plot(AX(2), binnedToPlot(:,1), stimsAfterLinearFilter);
        h.Color = COLOR_MAP(3, :);
        h.Color(4) = 0.45;  % 55% transparent
        h.LineStyle = '--';

        %G filter
        h = plot(AX(2), binnedToPlot(:,1), stimsAfterGenerator);
        h.Color = COLOR_MAP(2, :);
        h.Color(4) = 0.70;  % 30% transparent
        h.LineWidth = 1.5;

        legend('Raw stim values', ...
                'Aps', ...
                sprintf('PSTH (%.2f ms)', curBinSize/10), ...
                sprintf('After %s linear filter', ...
                    iif(Simulation.UsingSTA,'STA','STC')), ...
                'After Generator');
        
        %% fit measurement
        SSEk = curNeuron.Rate{iBinSize}.SSEk; %lower is less err
        SSEg = curNeuron.Rate{iBinSize}.SSEg; %lower is less err
                
        %the similarity between two signals, we only need zero lag
        Cork = abs(xcorr(psth,stimsAfterLinearFilter,0,'coeff')); % 1 if are equal
        Corg = abs(xcorr(psth,stimsAfterGenerator,0,'coeff')); % 1 if are equal

        yText = minValue-(minValue/10);
        
        
        ha = axes('Position',PositionOfSubplot,'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        
        text(0, 0.05, ...
            sprintf(['SSE after %s kernel: %.2f;\nSSE after Generator: %.2f;\n' ...
            'Cor.k: %.4f;\nCor.g: %.4f'], ...
                iif(Simulation.UsingSTA,'STA','STC'), ...
                SSEk, SSEg, Cork, Corg), ...
        'HorizontalAlignment' ,'left','VerticalAlignment', 'bottom');
        
    CreateTitleForSubplots(['\bf ' title]);
    
    r = 150; %pixels pre inch
    set(hf, 'PaperUnits', 'inches');
    set(hf, 'PaperPosition', [0 0 2880 1620]/r); %x_width=10cm y_width=15cm

    saveas(iBinSize, ['QnD_RvsRest_' MODE '_BinSize_' ...
        sprintf('%d', floor(curBinSize)) ...
        '_UsingSTA_' num2str(Simulation.UsingSTA)], 'png');
