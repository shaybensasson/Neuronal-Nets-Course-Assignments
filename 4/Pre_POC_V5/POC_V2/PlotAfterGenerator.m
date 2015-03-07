close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

load(['AfterGenerator_' MODE '.mat'])
    
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

NEURONS = length(Simulation.Neuron);

%iBinSize = 4;
%iNeuron = 2;

 for iBinSize=1:numel(Simulation.RATE_BIN_SIZES)
    curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    
    figure( 'Name', sprintf('After Generator (%s), Bin size: %.2f', ...
        MODE, curBinSize));
    
    for iNeuron=1:NEURONS
        subplot(2,2,iNeuron);
        
        curNeuron = Simulation.Neuron{iNeuron};

        times = curNeuron.Data(:,1);
        stimValues = curNeuron.Data(:,2);
        aps = curNeuron.Data(:,3);

        rateData = curNeuron.Rate{iBinSize}.Data;
        timesOfPSTH = rateData(:,1);
        psth = rateData(:,2);
        stimsAfterLinearFilter = rateData(:,3);
        stimsAfterGenerator = rateData(:, 4);

        %% plot

        map = [43, 87, 154; ... %blue: 1
            32, 162, 58; ... %green: 2
            0, 0, 0; ... % black: 3
            175, 185, 22; ... %yello for stims: 4
            201, 67, 67]; %red for aps: 5

        map = map./255;

        TIME_TO_PLOT = 10*Simulation.TICKS_IN_SECOND; %~ num of records
        START_FROM_INDEX = 1000;
        xStimFilter = times>=times(START_FROM_INDEX) & times<=times(START_FROM_INDEX)+TIME_TO_PLOT;
        timesToPlot = times(xStimFilter);
        dataToPlot = [times stimValues aps];
        dataToPlot = dataToPlot(xStimFilter, :);
        stimValues = dataToPlot(:,2);
        aps = dataToPlot(:,3);

        xRateFilter = rateData(:,1)>=times(START_FROM_INDEX) & ...
            rateData(:,1)<=times(START_FROM_INDEX)+TIME_TO_PLOT;
        binnedToPlot = rateData(xRateFilter, :);
        
        stackedToPlot = [rateData(:,1) psth stimsAfterLinearFilter stimsAfterGenerator];
        stackedToPlot = stackedToPlot(xRateFilter,:);

        %restore original times
        psth = stackedToPlot(:,2);
        stimsAfterLinearFilter = stackedToPlot(:,3);
        stimsAfterGenerator = stackedToPlot(:,4);

        %TODO: do we need this?
        %rate(isnan(rate(:, 2)), :) = [];


        %bind aps to top
        minValue = min(min(binnedToPlot(:,2:4)));
        maxValue = max(max(binnedToPlot(:,2:4)));
        
        aps = aps*maxValue;


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
        AX(2).YLim = [minValue maxValue];
        AX(2).YColor = map(1, :);

        %filter relevant data to plot
        AX(1).XLim = [timesToPlot(1) timesToPlot(end)];
        AX(2).XLim = [timesToPlot(1) timesToPlot(end)];


        %psth
        h = plot(AX(2), stackedToPlot(:,1), psth);
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
        
        %% fit measurement
        SSEk = curNeuron.Rate{iBinSize}.SSEk; %lower is less err
        SSEg = curNeuron.Rate{iBinSize}.SSEg; %lower is less err
                
        %the similarity between two signals, we only need zero lag
        Cork = xcorr(psth,stimsAfterLinearFilter,0,'coeff'); % 1 if are equal
        Corg = xcorr(psth,stimsAfterGenerator,0,'coeff'); % 1 if are equal

        yText = minValue-(minValue/10);
        text(timesToPlot(1), yText, ...
            sprintf(['SSE after STA kernel: %.2f;\nSSE after Generator: %.2f;\n' ...
            'Cor.k: %.4f;\nCor.g: %.4f'], ...
                SSEk, SSEg, Cork, Corg));
            
        
    end
    
 end
 
    





        

