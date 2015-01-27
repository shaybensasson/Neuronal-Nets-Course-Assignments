close all;

if (~exist('StimTimeRep','var'))
    load('fixedData.mat');
    %load('Data.mat');
end

%{
TT(i): contains data recorded from neuron i.
TT(i).sp: the times of the APs
a. times are in 1/10,000 second
b. from the onset of the experiment

*Rep: repeating stimuls every 200s in 30Hz (30 light-flickers per sec)  
->StimulusRep - stimulus 6000 light values that one of them appear each
    1/30 sec in the 200s repeating window
->StimTimeRep(t) = onset of the stimulus (prevWindowStart each repetition window),
    has to be normalized into 200 diffs, currently 100 diffs.

*NonRep: stimuls every 100s in 30Hz (30 light-flickers per sec)  
->StimulusNonRep - stimulus 3000 light values that one of them appear each
    1/30 sec in the 100s non repeating window
->StimTimeNonRep(t) = onset of the stimulus (prevWindowStart new stimulus window),
    alreay normalized into 100 diffs
->From StimTimeNonRep(1) till StimTimeNonRep(2), there were
    StimulusNonRep(1:3000) stimulus light values
    ->From StimTimeNonRep(3) till StimTimeNonRep(3), there were
        StimulusNonRep(3001:6000) stimulus light values
    ->and so prevWindowStart
%}

figureEx('Custom', 'Maximize');

%sampling rate is 1/10000, we multiply by 3 so it'd be 1/30000, because
%10000/30 is not an integer
TIME_FACTOR = 3; 
SAMPLING_TIME_RESOLUTION = 10000;
TICKS_IN_SECOND=TIME_FACTOR*SAMPLING_TIME_RESOLUTION; %ex. 30000 per sec
REPETITION_WINDOW_IN_SECS = 200; %in secs

NUM_OF_NEURONS=4; 
REPETITION_WINDOW=REPETITION_WINDOW_IN_SECS*TICKS_IN_SECOND; %200*3*10^4=6*10^6

%normalized into 200 diffs, currently 100 diffs
StimTimeRepFixed = StimTimeRep(1:2:end);
StimTimeRepFixed = StimTimeRepFixed*TIME_FACTOR;
%StimTimeRepFixed(2)-StimTimeRepFixed(1) ~= 200*10^4  *TIME_FACTOR

%how many repetitions do we have?
NUM_OF_REPETITIONS = length(StimTimeRepFixed);

STIMULI_VALUES = length(StimulusRep);

BIN_SIZE=100000;
NUM_OF_BINS =REPETITION_WINDOW/BIN_SIZE; %60 bins
bins = 0:BIN_SIZE:REPETITION_WINDOW;
indexesOfStimulusTime = 1:REPETITION_WINDOW;

maxHz = 0; %Max Hz for graphs
maxBinCount = 0; %Max bin count for graphs
handles = zeros(NUM_OF_NEURONS, 1); %used later to update graphs
%% plot neurons data
for iNeuron=1:NUM_OF_NEURONS
    
    subplot(NUM_OF_NEURONS,1,iNeuron);
    title(sprintf('Neuron #%d (Bin size/TicksPerSec = %d/%d=%0.2f secs)', ...
        iNeuron, BIN_SIZE,TICKS_IN_SECOND,BIN_SIZE/TICKS_IN_SECOND));
    
    hold on
    ylabel('Firing rate (Hz)')
    xlabel('Time of repetition window (s)')

    set(gca,'XTick',0:6*10^5:6*10^6);
    set(gca,'XTickLabel',0:20:REPETITION_WINDOW_IN_SECS);
    
    timeOfAPs=TTRep(iNeuron).sp*TIME_FACTOR;
    indexesPerNeuron = 1:length(timeOfAPs);
    accAPs = zeros(REPETITION_WINDOW, 1);
        
    for iIteration=1:NUM_OF_REPETITIONS
        windowStart = StimTimeRepFixed(iIteration);
        nextWindowStart = windowStart+REPETITION_WINDOW;
        
        %{
        Change from time scale into discrete indexes and 
        filter by AP occurences
        
        logical() is required to transform binary vector to a valid filter.
        %}
        filter = logical(timeOfAPs(:,1) >= windowStart & timeOfAPs(:,1) < nextWindowStart);
        indexesOfAPs = indexesPerNeuron(filter);
        
        
        %we have any APs
        if(~isempty(indexesOfAPs))
            %normalize the APs
            normalizedAPs = timeOfAPs(indexesOfAPs)-windowStart;
            
            curAPs = zeros(REPETITION_WINDOW, 1);
            curAPs(normalizedAPs) = 1;
            accAPs = accAPs + curAPs;
        end
    end %for ITERATIONS
    
    %plot historam for neuron
    filter = accAPs(:,1) > 0;
    spikeTimes = indexesOfStimulusTime(filter)';
    accAPsAtLeastOne = accAPs(filter);
    
    %{
    fprintf('[n=%d] sum accAPs=%d ...\n', ...
        iNeuron, sum(accAPs));
    %}
    
    [bincounts,binIndex] = histc(spikeTimes,bins);
    sumByBins = accumarray(binIndex,accAPsAtLeastOne, [length(bins) 1]);
    bar(bins,sumByBins,'histc');
    curMaxBinCount = max(sumByBins);
    maxBinCount = max(curMaxBinCount, maxBinCount);
    
    curHz = curMaxBinCount/(BIN_SIZE/TICKS_IN_SECOND);
    
    maxHz = max(curHz, maxHz);
    handles(iNeuron)=gca;
end


for iNeuron=1:NUM_OF_NEURONS
    h = handles(iNeuron);
    axis(h, [0 REPETITION_WINDOW 0 maxBinCount]);
    %ylim([0 maxBinCount]);
    med = round(maxBinCount/2);
    set(h,'YTick',[0 med maxBinCount]);
    med = round(maxHz/2);
    %set(h,'YTickLabel',{' ';'1';'2'});
    set(h,'YTickLabel',{' ';num2str(med);num2str(floor(maxHz))});
end
    
        
beep('on');
