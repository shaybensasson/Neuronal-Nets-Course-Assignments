if (~exist('TTNonRep','var'))
    load('FixedData.mat')
end

SECONDS_IN_WINDOW = 100;
TICKS_IN_SECOND = 10000;
TICKS_IN_WINDOW = TICKS_IN_SECOND * SECONDS_IN_WINDOW;
STIMULI_PER_SECOND = 1/30;
STIMULUS_EACH_TICKS = floor(TICKS_IN_SECOND*STIMULI_PER_SECOND); %333
STIMULI_PER_WINDOW = SECONDS_IN_WINDOW / STIMULI_PER_SECOND; %3000 stims for 100 secs

NEURONS = length(TTNonRep);
ITERATIONS = length(StimTimeNonRep);
ITERATIONS = ITERATIONS-2; %we ignore the last 2 chunks because it has partial stimuli values

SAFETY_SIZE_SUFFIX = 1000; %actually the max is 3725

        
Neuron = cell(1,length(TTNonRep));
for iNeuron = 1:length(TTNonRep)
    

    simulation.all = NaN((STIMULI_PER_WINDOW+SAFETY_SIZE_SUFFIX)*ITERATIONS, 3);
    lastIterationIndex = 0;
    lastIterationSimuliIndex = 0;
    
    for iIteration = 2:ITERATIONS
        fprintf('[N:#%i] processing iteration #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 

        StimTimeNonRepNorm=StimTimeNonRep-StimTimeNonRep(iStart);
        iteration.time= (StimTimeNonRepNorm(iStart)+1:1:StimTimeNonRepNorm(iEnd)+1)'; %first 100 secs

        actualWindowLength = StimTimeNonRep(iEnd)-StimTimeNonRep(iStart)+1;

        
        %window stim values by time, warming it up, to improve prefs
        stimTimes = 1:STIMULUS_EACH_TICKS:actualWindowLength;
        nextIterationSimuliIndex = lastIterationSimuliIndex + length(stimTimes) + 1;
        globalStimValues = StimulusNonRep(lastIterationSimuliIndex+1:nextIterationSimuliIndex-1,:);
        lastIterationSimuliIndex = nextIterationSimuliIndex+1;
        
        stimValues = globalStimValues(:);
        iteration.stimuli = NaN(actualWindowLength,1);
        iteration.stimuli(stimTimes)=stimValues;

        timeOfAPs=TTNonRep(iNeuron).sp; %time of Aps

        onset = StimTimeNonRep(iStart)+1;
        next_onset = StimTimeNonRep(iEnd)+1;
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





