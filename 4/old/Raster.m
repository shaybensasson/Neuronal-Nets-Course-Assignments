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
->StimTimeRep(t) = onset of the stimulus (on each repetition window),
    has to be normalized into 200 diffs, currently 100 diffs.

*NonRep: stimuls every 100s in 30Hz (30 light-flickers per sec)  
->StimulusNonRep - stimulus 3000 light values that one of them appear each
    1/30 sec in the 100s non repeating window
->StimTimeNonRep(t) = onset of the stimulus (on new stimulus window),
    alreay normalized into 100 diffs
->From StimTimeNonRep(1) till StimTimeNonRep(2), there were
    StimulusNonRep(1:3000) stimulus light values
    ->From StimTimeNonRep(3) till StimTimeNonRep(3), there were
        StimulusNonRep(3001:6000) stimulus light values
    ->and so on
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
ITERATIONS = length(StimTimeRepFixed);

for iNeuron=1:NUM_OF_NEURONS
    subplot(NUM_OF_NEURONS,1,iNeuron);
    hold on;
    title(sprintf('Neuron #%d', iNeuron));
    ylabel('# of Iteration')
    xlabel('Time of repetition window (s)')
    ylim([0 ITERATIONS]);
    set(gca,'XTickLabel',0:20:REPETITION_WINDOW_IN_SECS);
    
    
    
    for iIteration=1:ITERATIONS
        timeOfAPs = TTRep(1,iNeuron).sp * TIME_FACTOR;
        
        %current window/rep of 200 secs, starting time
        windowStart = StimTimeRepFixed(iIteration);
                
        %look until next window
        nextWindowStart= windowStart+REPETITION_WINDOW;
        
        indexes = 1:length(timeOfAPs);
        ind = indexes(timeOfAPs(:) >= windowStart & timeOfAPs(:) < nextWindowStart);
        
        %get all APs in the window
        normalizedAPs = timeOfAPs(ind)-windowStart;
        %get back to 10000 scale
        normalizedAPs = normalizedAPs / TIME_FACTOR;
              
        scatter(normalizedAPs, ...
            ones(length(normalizedAPs),1)*iIteration,1,'b');
    end 
end