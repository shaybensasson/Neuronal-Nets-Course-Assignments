close all;

PSTH_BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq
             Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec          
             Simulation.STIMULUS_EACH_TICKS*5; ... ~150 msec
             Simulation.STIMULUS_EACH_TICKS*10; ... ~330 msec
             Simulation.STIMULUS_EACH_TICKS*15; ... ~500 msec
             Simulation.STIMULUS_EACH_TICKS*30]; %1000 msec = 1 sec  
curBinSize = PSTH_BIN_SIZES(2);

data = Simulation.Neuron{2}.Data;

times = data(:,1);
stimValues = data(:,2);
aps = data(:,3);

rateData = [times aps];

%bins = onset:curBinSize:next_onset-1;
bins = times(1):curBinSize:times(end);
[~,binIndex] = histc(times,bins);

%insert any unbinned data to last bin
binIndex(binIndex==0)=max(binIndex);

%group by times and sum aps
grp_psth = accumarray(binIndex, aps, [length(bins) 1], @nansum, NaN);

binned = [ceil(bins')-times(1)+1 grp_psth bins'];

firstApTime = data(~isnan(aps) & aps>0, 1);
firstApTime = firstApTime(1);

%throw bins before first ap
binned(binned(:,2) == 0 & binned(:,3)<firstApTime, :) = [];

%throw last bin, it has garbage unbinned data (it's acc is NaN)
binned = binned(1:end-1,:); 
psth = binned(:, 2);

NORMALIZE = 1;
NORMALIZE_BY_MAX = 1;

psth = normalize(psth, NORMALIZE, NORMALIZE_BY_MAX);
stimValues = normalize(stimValues, NORMALIZE, NORMALIZE_BY_MAX);

%% plot
figure;

hold on;
[AX,H1,H2] = plotyy(times, stimValues, times, aps, 'plot', 'plot');
%hold(AX(1));
set(AX,'NextPlot','add')

H1.LineStyle = 'none'; H1.Marker = 'o';
AX(1).YLim = [min(stimValues) max(stimValues)];
H2.LineStyle = 'none'; H2.Marker = 'o';
AX(2).YLim = [min(psth) max(psth)];

%filter relevant data to plot
TIME_TO_PLOT = 30*TICKS_IN_SECOND;
timesToPlot = times(times<=times(1)+TIME_TO_PLOT);
AX(1).XLim = [timesToPlot(1) timesToPlot(end)];
AX(2).XLim = [timesToPlot(1) timesToPlot(end)];

rate = [binned(:,3) psth];
rate(isnan(psth(:)), :) = [];

plot(AX(2), rate(:,1), rate(:,2));


