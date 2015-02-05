close all;

ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)) || ...
        Simulation.Phase < CONSTANTS.PHASES.GENERATOR)
    clearvars -except MODE;
    load(['AfterGenerator_' MODE '.mat'])
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

NEURONS = length(Simulation.Neuron);
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;
SECONDS_OF_RATE_TO_DISPLAY = Simulation.SECONDS_OF_RATE_TO_DISPLAY;
ITERATIONS = Simulation.ITERATIONS;
PSTH_BIN_SIZES = Simulation.PSTH_BIN_SIZES;
             
STIMULI_PER_WINDOW = Simulation.STIMULI_PER_WINDOW;
%store for later usage


%for iBinSize=4:6
for iBinSize=1:numel(PSTH_BIN_SIZES)
    curBinSize = PSTH_BIN_SIZES(iBinSize);
    
    
    figure( 'Name', sprintf('RvsRest (%s), Bin size: %.2f', ...
        MODE, curBinSize));

    %for iNeuron=3:3
    for iNeuron=1:NEURONS
        %figure( 'Name', sprintf('RvsRest (%s) of N#%d, Bin size: %.2f', ...
         %   MODE, iNeuron, curBinSize));
    
        subplot(2,2,iNeuron);
        curNeuron = Simulation.Neuron{iNeuron};
        fprintf('[Bsz,N:#%d,#%d] ...\n', iBinSize, iNeuron);
        
        data = curNeuron.PSTH{iBinSize}.RVsRest;
        times = data(:,1);
        PSTH = data(:,2);
        stimsAfterLinearFilter = data(:,3);
        stimsAfterGenerator = data(:,4);
        
        hold on;
           
        h=plot(times,stimsAfterLinearFilter,'k');
        h.Color(4) = 0.3;  % 70% transparent
        %h = plot(times,stimValues, 'b');
        h = plot(times,PSTH, 'b');
        h.Color(4) = 0.5;  % 50% transparent
        h = plot(times,stimsAfterGenerator,'g');
        h.Color(4) = 0.8;  % 20% transparent

        title(sprintf('Neuron #%d', iNeuron));
        legend( ...
            'After STA Linear Filter', ...
            sprintf('PSTH (Bin Size = %.2f ms)', ...
                curBinSize/Simulation.TICKS_IN_SECOND*1000), ...
            'After Generator (Rest)');
        xlim([0 times(end)]);

        %{
        set(gca,'XTickLabel',sprintf('%1.0f|',...
            0:SECONDS_IN_WINDOW:times(end)/TICKS_IN_SECOND));
        %}

        xlabel('Time (s)');
        ylabel('Normalized Values');

        hold off;

        SSEk = curNeuron.PSTH{iBinSize}.SSEk; %lower is less err
        SSEg = curNeuron.PSTH{iBinSize}.SSEg; %lower is less err
                
        %the similarity between two signals, we only need zero lag
        %CC = xcorr(PSTH,stimsAfterGenerator,0,'coeff'); % 1 if are equal

        text(0, 0.9, ...
            sprintf('SSE after STA kernel: %.2f;\nSSE after Generator: %.2f;', ...
                SSEk, SSEg));

    end %for iNeuron
end %for iBinSize
        
beep('on');

