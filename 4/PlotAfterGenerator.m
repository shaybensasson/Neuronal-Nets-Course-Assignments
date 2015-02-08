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

timesOfPSTH = binned(:,4);
psth = binned(:,2);
stimsAfterLinearFilter = binned(:,3);
%ORIGINAL_stimsAfterGenerator = stimsAfterGenerator; %for multiple tests
stimsAfterGenerator = ORIGINAL_stimsAfterGenerator;

%% plot

map = [43, 87, 154; ... %blue: 1
    32, 162, 58; ... %green: 2
    0, 0, 0; ... % black: 3
    175, 185, 22; ... %yello for stims: 4
    201, 67, 67]; %red for aps: 5

map = map./255;

hf = figure;
TIME_TO_PLOT = 10*TICKS_IN_SECOND; %~ num of records
START_FROM_INDEX = 1000;
xStimFilter = times>=times(START_FROM_INDEX) & times<=times(START_FROM_INDEX)+TIME_TO_PLOT;
timesToPlot = times(xStimFilter);
dataToPlot = [times stimValues aps];
dataToPlot = dataToPlot(xStimFilter, :);
stimValues = dataToPlot(:,2);
aps = dataToPlot(:,3);

xRateFilter = binned(:,4)>=times(START_FROM_INDEX) & ...
    binned(:,4)<=times(START_FROM_INDEX)+TIME_TO_PLOT;
binnedToPlot = binned(xRateFilter, :);

%restore original times
stackedToPlot = [binned(:,4) psth stimsAfterLinearFilter stimsAfterGenerator];
stackedToPlot = stackedToPlot(xRateFilter,:);

psth = stackedToPlot(:,2);
stimsAfterLinearFilter = stackedToPlot(:,3);
stimsAfterGenerator = stackedToPlot(:,4);

%TODO: do we need this?
%rate(isnan(rate(:, 2)), :) = [];


%bind aps to top
aps = aps*max(binnedToPlot(:,2));


hold on;
[AX,H1,H2] = plotyy(timesToPlot, stimValues, timesToPlot, aps, 'plot', 'plot');
%hold(AX(1));
set(AX,'NextPlot','add')
H1.LineStyle = 'none'; H1.Marker = '.'; 
H1.Color = map(4, :);  
H1.Color(4) = 0.30;  % 70% transparent

AX(1).YLim = [min(stimValues) max(stimValues)];
AX(1).YColor = map(4, :);

H2.LineStyle = 'none'; H2.Marker = 'o';
H2.Color = map(5, :);
H2.Color(4) = 0.30;  % 70% transparent
AX(2).YLim = [min(min(stackedToPlot(:,2:4))) max(max(stackedToPlot(:,2:4)))];
AX(2).YColor = map(1, :);

%filter relevant data to plot
AX(1).XLim = [timesToPlot(1) timesToPlot(end)];
AX(2).XLim = [timesToPlot(1) timesToPlot(end)];


%psth
h = plot(AX(2), stackedToPlot(:,1), stackedToPlot(:,2));
h.Color = map(1, :);
h.Color(4) = 0.60; % 30% transparent 
h.LineWidth = 1.5;
%K filter
h = plot(AX(2), stackedToPlot(:,1), stimsAfterLinearFilter);
h.Color = map(3, :);
h.Color(4) = 0.35;  % 65% transparent
h.LineStyle = '--';

%G filter
h = plot(AX(2), stackedToPlot(:,1), stimsAfterGenerator);
h.Color = map(2, :);
h.Color(4) = 0.70;  % 30% transparent
h.LineWidth = 1.5;

legend('Raw stim values', ...
        'Aps', ...
        sprintf('PSTH (%.2f ms)', curBinSize/10), ...
        'After linear filter', ...
        'After Generator');
        

