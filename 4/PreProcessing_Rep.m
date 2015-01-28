if (~exist('TTRep','var'))
    load('FixedData.mat')
    
    %fix diffs to be 200s instead of 100s, each time we load raw data
    StimTimeRep = StimTimeRep(1:2:end); 
end

SECONDS_IN_WINDOW = 200;
TICKS_IN_SECOND = 10000;
TICKS_IN_WINDOW = TICKS_IN_SECOND * SECONDS_IN_WINDOW;
STIMULI_PER_SECOND = 1/30;
STIMULUS_EACH_TICKS = floor(TICKS_IN_SECOND*STIMULI_PER_SECOND); %333
STIMULI_PER_WINDOW = SECONDS_IN_WINDOW / STIMULI_PER_SECOND; %6000 stims for 200 secs

NEURONS = length(TTRep);
ITERATIONS = length(StimTimeRep);

SAFETY_SIZE_SUFFIX = 10; %actually the max is 6008

%window stim values by time, warming it up, to improve prefs
globalStimValues = StimulusRep(1:STIMULI_PER_WINDOW);
globalStimValues = repmat(globalStimValues,2,1);

maxstimTimes = 0;

Neuron = cell(1,length(TTRep));
for iNeuron = 1:length(TTRep)
    

    simulation.all = NaN((STIMULI_PER_WINDOW+SAFETY_SIZE_SUFFIX)*ITERATIONS, 3);
    lastIterationIndex = 0;
    for iIteration = 2:ITERATIONS
        %fprintf('[N:#%i] processing iteration #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 

        StimTimeRepNorm=StimTimeRep-StimTimeRep(iStart);
        iteration.time= (StimTimeRepNorm(iStart)+1:1:StimTimeRepNorm(iEnd)+1)'; %first 200 secs

        actualWindowLength = StimTimeRep(iEnd)-StimTimeRep(iStart)+1;

        
        %copy global var
        stimValues = globalStimValues(:);
        iteration.stimuli = NaN(actualWindowLength,1);
        stimTimes = 1:STIMULUS_EACH_TICKS:actualWindowLength;
        maxstimTimes = max(maxstimTimes, length(stimTimes));
        iteration.stimuli(stimTimes)=stimValues(1:length(stimTimes));

        timeOfAPs=TTRep(iNeuron).sp; %time of Aps

        onset = StimTimeRep(iStart)+1;
        next_onset = StimTimeRep(iEnd)+1;
        %{
                Change from time scale into discrete indexes and 
                filter by AP occurences

                logical() is required to transform binary vector to a valid filter.
                %}
        indexes = 1:length(timeOfAPs);
        filter = logical(timeOfAPs(:) >= onset & timeOfAPs(:) < next_onset);
        indexesOfAPs = indexes(filter);

        iteration.APs = zeros(actualWindowLength,1);
        %normalize the APs, to the start of window
        APsInWindow = timeOfAPs(indexesOfAPs)-onset + 1; %the index is 1 based
        iteration.APs(APsInWindow)=1;

        %create fixed data set
        res = NaN(actualWindowLength, 3);
        res(:,1) = iteration.time + onset;
        res(:,2) = iteration.stimuli;
        res(:,3) = iteration.APs;

        iteration.all = res(res(:,2) > 0 | res(:,3) > 0,:);
        nextIterationIndex = lastIterationIndex + length(iteration.all) + 1;
        simulation.all(lastIterationIndex+1:nextIterationIndex-1,:) = iteration.all;
        lastIterationIndex = nextIterationIndex+1;
    end %iIteration

    %remove NaN rows
    Neuron{iNeuron} = simulation.all(any(~isnan(simulation.all),2),:);
end %iNeuron