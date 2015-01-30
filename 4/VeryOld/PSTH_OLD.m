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

%normalized into 200 diffs, currently 100 diffs
StimTimeRepEvery200 = StimTimeRep(1:2:end);


%sampling rate is 1/10000, we multiply by 3 so it'd be 1/30000, because
%10000/30 is not an integer
TIME_FACTOR = 3; 
SAMPLING_TIME_RESOLUTION = 10000;
TICKS_IN_SECOND=TIME_FACTOR*SAMPLING_TIME_RESOLUTION; %ex. 30000 per sec
REPETITION_WINDOW_IN_SECS = 200; %in secs

NUM_OF_NEURONS=4; 
REPETITION_WINDOW=REPETITION_WINDOW_IN_SECS*TICKS_IN_SECOND; %200*3*10^4=6*10^6

%normalized into 200 diffs, currently 100 diffs
StimTimeRepFixed = StimTimeRep(1:2:end)*TIME_FACTOR;
%StimTimeRepFixed(2)-StimTimeRepFixed(1) ~= 200*10^4  *TIME_FACTOR

%how many repetitions do we have?
ITERATIONS = length(StimTimeRepFixed);

%convert to 1/30000 sec
noOfstimuliPerRun = 6000;

BIN_SIZE=100000;

numOfBins =REPETITION_WINDOW/BIN_SIZE; %60 bins

PLOT_RATE = 1;


StimTimeRepFixed = StimTimeRep(1:2:end)*TIME_FACTOR;

for iNeuron=1:NUM_OF_NEURONS
    subplot(NUM_OF_NEURONS+1,1,iNeuron);
    hold on;
    title(sprintf('Neuron #%d', iNeuron));
    
    psth = zeros(numOfBins,1);
    for iIteration=1:ITERATIONS
        timeOfAPs = TTRep(1,iNeuron).sp * TIME_FACTOR;
        
        
        windowStart = StimTimeRepFixed(iIteration);
        
        nextWindowStart= windowStart+REPETITION_WINDOW;
        
        
        indexes = 1:length(timeOfAPs);
        ind = indexes(timeOfAPs(:) >= windowStart & timeOfAPs(:) < nextWindowStart);
        
        
                
        
        fromBin=0;
        %get all APs in the window
        normalizedAPs= timeOfAPs(ind)-windowStart;
        for iBin=1:numOfBins
            toBin = fromBin+BIN_SIZE;
            indexes = 1:length(normalizedAPs);
            indBin = indexes(normalizedAPs(:) >= fromBin & normalizedAPs(:) < toBin);
            psth(iBin) = psth(iBin)+length(normalizedAPs(indBin));
            fromBin = toBin;
        end
    end
    
    if(PLOT_RATE==1)
        bar(0:BIN_SIZE:(REPETITION_WINDOW-BIN_SIZE),(psth/(noOfstimuliPerRun*length(StimTimeRepFixed)*BIN_SIZE))*REPETITION_WINDOW)
        
        ylabel('Rate(Hz)');
        axis([-1 REPETITION_WINDOW 0 max((psth/(noOfstimuliPerRun*length(StimTimeRepFixed)*BIN_SIZE))*REPETITION_WINDOW)])
    else
        bar(0:BIN_SIZE:(REPETITION_WINDOW-BIN_SIZE),psth)
        
        axis([-1 REPETITION_WINDOW 0 max(psth)])
    end
    xlabel(['Time (12000000= 400 sec) , Bin size in sec: ' num2str(BIN_SIZE/TICKS_IN_SECOND) ] )
    
    
end