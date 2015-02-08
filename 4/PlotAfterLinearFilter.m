close all;

PSTH_BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq (1)
             Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec (2)          
             Simulation.STIMULUS_EACH_TICKS*5; ... ~150 msec (3)
             Simulation.STIMULUS_EACH_TICKS*10; ... ~330 msec (4)
             Simulation.STIMULUS_EACH_TICKS*15; ... ~500 msec (5)
             Simulation.STIMULUS_EACH_TICKS*30]; %1000 msec = 1 sec (6)
curBinSize = PSTH_BIN_SIZES(4);

data = Simulation.Neuron{2}.Data;

times = data(:,1);
stimValues = data(:,2);
aps = data(:,3);
stimValuesAfterK = data(:,4);

rateData = [times aps];

%bins = onset:curBinSize:next_onset-1;
bins = times(1):curBinSize:times(end);
[~,binIndex] = histc(times,bins);

%insert any unbinned data to last bin
binIndex(binIndex==0)=max(binIndex);

%group by times and sum aps
grp_psth = accumarray(binIndex, aps, [length(bins) 1], @nansum, NaN);
%group by times and mean stim vals
meanIgnoreNaNs = @(vector) mean(vector(~isnan(vector(:))));
grp_afterK = accumarray(binIndex, stimsAfterLinearFilter, [length(bins) 1], meanIgnoreNaNs, NaN);

binned = [ceil(bins')-times(1)+1 grp_psth grp_afterK bins'];

firstApTime = data(~isnan(aps) & aps>0, 1);
firstApTime = firstApTime(1);

%remove NaNs
binned(isnan(binned(:, 2)) | isnan(binned(:, 3)), :) = [];

%throw bins before first ap
binned(binned(:,2) == 0 & binned(:,4)<firstApTime, :) = [];

%throw last bin, it has garbage unbinned data (it's acc is NaN)
binned = binned(1:end-1,:); 
psth = binned(:, 2);
meanStimAfterK = binned(:, 3);

NORMALIZE = 1;
NORMALIZE_BY_MAX = 1;

%TODO: normalize while storing in curNeuron.Data
psth = normalize(psth, NORMALIZE, NORMALIZE_BY_MAX);
binned(:, 2) = psth;
meanStimAfterK = normalize(meanStimAfterK, NORMALIZE, NORMALIZE_BY_MAX);
binned(:, 3) = meanStimAfterK;

stimValues = normalize(stimValues, NORMALIZE, NORMALIZE_BY_MAX);
data(:,2) = stimValues;
stimValuesAfterK = normalize(stimValuesAfterK, NORMALIZE, NORMALIZE_BY_MAX);
data(:,4) = stimValuesAfterK;

%% plot
hf = figure;
TIME_TO_PLOT = 10*TICKS_IN_SECOND; %~ num of records
START_FROM_INDEX = 1000;
xStimFilter = times>=times(START_FROM_INDEX) & times<=times(START_FROM_INDEX)+TIME_TO_PLOT;
timesToPlot = times(xStimFilter);
dataToPlot = [times stimValues aps stimValuesAfterK];
dataToPlot = dataToPlot(xStimFilter, :);
stimValues = dataToPlot(:,2);
aps = dataToPlot(:,3);
stimValuesAfterK = dataToPlot(:,4);

xRateFilter = binned(:,4)>=times(START_FROM_INDEX) & ...
    binned(:,4)<=times(START_FROM_INDEX)+TIME_TO_PLOT;
binnedToPlot = binned(xRateFilter, :);

%restore original times
rate = [binned(:,4) psth];
rate = rate(xRateFilter,:);

%TODO: do we need this?
%rate(isnan(rate(:, 2)), :) = [];


%bind aps to top
aps = aps*max(binnedToPlot(:,2));


hold on;
[AX,H1,H2] = plotyy(timesToPlot, stimValues, timesToPlot, aps, 'plot', 'plot');
%hold(AX(1));
set(AX,'NextPlot','add')
H1.LineStyle = 'none'; H1.Marker = '.'; 
H1.Color(4) = 0.30;  % 65% transparent

AX(1).YLim = [min(min(stimValues, stimValuesAfterK)) max(max(stimValues, stimValuesAfterK))];
H2.LineStyle = 'none'; H2.Marker = 'o';
AX(2).YLim = [min(binnedToPlot(:,2)) max(binnedToPlot(:,2))];

%filter relevant data to plot
AX(1).XLim = [timesToPlot(1) timesToPlot(end)];
AX(2).XLim = [timesToPlot(1) timesToPlot(end)];


plot(AX(2), rate(:,1), rate(:,2), 'b');

h = plot(AX(1), timesToPlot, stimValuesAfterK, 'k');
h.Color(4) = 0.35;  % 65% transparent

legend('Raw stim values', ...
        'After linear filter', ...
        'Aps', ...
        sprintf('PSTH (%.2f ms)', curBinSize/10));

