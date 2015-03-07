if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,'NonRep')))
    clear all;
    load('PreProcessed_NonRep.mat')
end

ConstantsHeader();
SAVE_MAT_FILE = 1;

NEURONS = length(Simulation.Neuron);
ITERATIONS = Simulation.ITERATIONS;
StimTimeNonRep = Simulation.StimTimeNonRep;
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;

%BIN_SIZE=10/3*TICKS_IN_SECOND; %3.3333*10^4
%BIN_SIZE=1000; %0.1 sec
%BIN_SIZE=Simulation.STIMULUS_EACH_TICKS; %Minimal BinSize = sampling freq
BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq
            Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec
            0.5*Simulation.TICKS_IN_SECOND; ... 0.5 sec
            Simulation.TICKS_IN_SECOND]; %1 sec

BIN_SIZE_TO_PLOT = BIN_SIZES(1);
indexesOfTime = 1:TICKS_IN_WINDOW; %used for x axis values

maxHz = 0; %Max Hz for graphs
maxBinCount = 0; %Max bin count for graphs
handles = zeros(NEURONS, 1); %used later to update graphs

Simulation.Phase = CONSTANTS.PHASES.PSTH;
Simulation.PSTH_BIN_SIZES = BIN_SIZES;

%% plot neurons data
for iNeuron=1:NEURONS
    fprintf('[N:#%i] ...\n', iNeuron);
    
    accAPs = zeros(TICKS_IN_WINDOW, 1);
    
    subplot(NEURONS,1,iNeuron);
    title(sprintf('Neuron #%d (Bin size/TicksPerSec = %.2f/%d=%.2f secs)', ...
        iNeuron, BIN_SIZE_TO_PLOT,TICKS_IN_SECOND,BIN_SIZE_TO_PLOT/TICKS_IN_SECOND));
    
    hold on
    
    ylabel('Firing rate (Hz)')
    xlabel('Time of repetition window (s)')

    set(gca,'XTick',0:2*10^5:2*10^6);
    set(gca,'XTickLabel',0:20:SECONDS_IN_WINDOW);
    
    for iIteration=2:ITERATIONS
        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 

        times = Simulation.Neuron{iNeuron}.Data(:,1);
        
        onset = StimTimeNonRep(iStart);
        next_onset = StimTimeNonRep(iEnd);
        filter = logical(times(:) >= onset & times(:) < next_onset);
        windowData = Simulation.Neuron{iNeuron}.Data(filter, :);
                
        delta = zeros(length(windowData), 3);
        delta(:,1) = onset;
        windowData=windowData-delta; %normalize
        
        %throw extra
        filter = logical(windowData(:,1) <= TICKS_IN_WINDOW);
        windowData = windowData(filter, :);
        
        %only APs
        filter = logical(windowData(:,3) == 1);
        windowData = windowData(filter, :);
        
        %we have any APs
        if(~isempty(windowData))
            curAPs = zeros(TICKS_IN_WINDOW, 1);
            curAPs(windowData(:,1)) = 1;
            accAPs = accAPs + curAPs;
        end
        
    end %for ITERATIONS
    
    maxy = ITERATIONS;
    axis([0 TICKS_IN_WINDOW 0 maxy]);
    med = round(maxy/2);
    set(gca,'YTick',[0 med maxy]);
    
    %plot histogram for neuron
    spikeTimes = indexesOfTime';
    
    %{
    fprintf('[n=%d] sum accAPs=%d ...\n', ...
        iNeuron, sum(accAPs));
    %}
    
    for iBinSize=1:numel(BIN_SIZES)

        binSize=BIN_SIZES(iBinSize);
        bins = 0:binSize:TICKS_IN_WINDOW;
    
        [bincounts,binIndex] = histc(spikeTimes,bins);
        sumByBins = accumarray(binIndex,accAPs, [length(bins) 1]);
        Simulation.Neuron{iNeuron}.PSTH{iBinSize} = [bins' sumByBins];

        if (binSize==BIN_SIZE_TO_PLOT) 
            bar(bins,sumByBins,'histc');
            curMaxBinCount = max(sumByBins);
            maxBinCount = max(curMaxBinCount, maxBinCount);
    
            curHz = curMaxBinCount/(BIN_SIZE_TO_PLOT/TICKS_IN_SECOND);
    
            maxHz = max(curHz, maxHz);
            handles(iNeuron)=gca;
        end
    end
    
end %for iNeuron

for iNeuron=1:NEURONS
    h = handles(iNeuron);
    axis(h, [0 TICKS_IN_WINDOW 0 maxBinCount]);
    ylim([0 maxBinCount]);
    med = round(maxBinCount/2);
    set(h,'YTick',[0 med maxBinCount]);
    med = round(maxHz/2);
    set(h,'YTickLabel',{' ';'1';'2'});
    set(h,'YTickLabel',{' ';num2str(med);num2str(floor(maxHz))});
end

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save('AfterPSTH_NonRep.mat', 'Simulation');
end
        
beep('on');
