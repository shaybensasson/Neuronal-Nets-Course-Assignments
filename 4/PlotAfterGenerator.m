close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();


%choose Rep or NonRep
MODE = 'NonRep';

load(['MatFiles\AfterGenerator_' MODE '.mat'])


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

%duration in seconds to display when ploting multiple rates
LENGTH_OF_PLOT_IN_SECS = 10;
LENGTH_OF_PLOT_IN_TICKS = LENGTH_OF_PLOT_IN_SECS*Simulation.TICKS_IN_SECOND; %~ num of records
START_PLOT_FROM_INDEX = 1;

COLOR_MAP = [43, 87, 154; ... %blue: 1
            32, 162, 58; ... %green: 2
            0, 0, 0; ... % black: 3
            175, 185, 22; ... %yello for stims: 4
            201, 67, 67]; %red for aps: 5

COLOR_MAP = COLOR_MAP./255;


%for iBinSize=1:1
for iBinSize=1:numel(Simulation.RATE_BIN_SIZES)
    curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    
    title = sprintf('R vs R_{Est} (%s), Bin size: %.2f, Using %s Filter', ...
        MODE, curBinSize, iif(Simulation.UsingSTA,'STA','STC'));
    
    hf = figure(iBinSize);
    hf.Name = title;
    
    %for iNeuron=2:2
    for iNeuron=1:NEURONS
        hs = subplot(2,2,iNeuron);
        
        PositionOfSubplot = hs.Position;
        
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
        
        
    end %iNeuron
    
    CreateTitleForSubplots(['\bf ' title]);
    
    r = 150; %pixels pre inch
    set(hf, 'PaperUnits', 'inches');
    set(hf, 'PaperPosition', [0 0 2880 1620]/r); %x_width=10cm y_width=15cm

    saveas(hf, ['RvsRest_' MODE '_BinSize_' ...
        sprintf('%d', floor(curBinSize)) ...
        '_UsingSTA_' num2str(Simulation.UsingSTA)], 'png');
    
 end %iBinSize
 
    





        

