clear all;
if (~exist('TTNonRep','var'))
    load('FixedData.mat')
end

ConstantsHeader();

SAVE_MAT_FILE = 1;

%secs in stimulus non-rep window
SECONDS_IN_WINDOW = 100;
TICKS_IN_SECOND = 10000;
TICKS_IN_WINDOW = TICKS_IN_SECOND * SECONDS_IN_WINDOW;
STIMULI_PER_SECOND = 1/30;
%how many ticks there are between stimuli.
STIMULUS_EACH_TICKS = TICKS_IN_SECOND*STIMULI_PER_SECOND; %333.333
STIMULI_PER_WINDOW = SECONDS_IN_WINDOW / STIMULI_PER_SECOND; %3000 stims for 100 secs

NEURONS = length(TTNonRep);
ITERATIONS = length(StimTimeNonRep); %NOTE: each ITERATION is different bulk (non-repeating window)
ITERATIONS = ITERATIONS-2; %we ignore the last 2 chunks because it has partial stimuli values

clear Simulation;
Simulation.Mode = CONSTANTS.MODES.NONREP;
Simulation.Phase = CONSTANTS.PHASES.PREPROCESSING;
Simulation.Neuron = cell(1,length(TTNonRep));
Simulation.StimTimeNonRep = StimTimeNonRep;
Simulation.ITERATIONS = ITERATIONS;
Simulation.SECONDS_IN_WINDOW = SECONDS_IN_WINDOW;
Simulation.TICKS_IN_WINDOW = TICKS_IN_WINDOW;
Simulation.TICKS_IN_SECOND = TICKS_IN_SECOND;
Simulation.STIMULI_PER_SECOND = STIMULI_PER_SECOND;
Simulation.STIMULUS_EACH_TICKS = STIMULUS_EACH_TICKS;
Simulation.STIMULI_PER_WINDOW = STIMULI_PER_WINDOW;

        
for iNeuron = 1:length(TTNonRep)
    
    fprintf('[N:#%i] ...\n', iNeuron);
    
    %allocate an array to store prepocessed data for neuron
    neuronData.all = NaN(STIMULI_PER_WINDOW*ITERATIONS, 3);
    lastIterationIndex = 0;
    lastIterationSimuliIndex = 0;
    
    for iIteration = 2:ITERATIONS
        %fprintf('[N:#%i] processing iteration #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1;
        
        %actual ticks in window limited by TICKS_IN_WINDOW
        dataTicks = TICKS_IN_WINDOW;

        %Normalize stim time to start from 1
        StimTimeNonRepNorm=StimTimeNonRep-StimTimeNonRep(iStart)+1;
        iteration.time= (StimTimeNonRepNorm(iStart):StimTimeNonRepNorm(iEnd)-1)'; %first 100 secs
        dataTicks = min(dataTicks,length(iteration.time)); 
        iteration.time = iteration.time(1:dataTicks); %trim
        
        %% get stim values of the current non-rep window        
        iteration.stimuli = NaN(dataTicks,1);
	
        %create stim times for every stim
        stimTimes = 1:floor(STIMULUS_EACH_TICKS):dataTicks;
        stimTimes = stimTimes(1:STIMULI_PER_WINDOW);
        %get stim values on non-rep window, using stim times
        nextIterationSimuliIndex = lastIterationSimuliIndex + length(stimTimes) + 1;
        stimValues = StimulusNonRep(lastIterationSimuliIndex+1:nextIterationSimuliIndex-1,:);
        lastIterationSimuliIndex = nextIterationSimuliIndex-1;
                
        iteration.stimuli(stimTimes)=stimValues;

        %% get aps of the current non-rep window        
        timeOfAPs=TTNonRep(iNeuron).sp; %time of Aps

        onset = StimTimeNonRep(iStart);
        next_onset = StimTimeNonRep(iEnd);
        %{
                Change from time scale into discrete indexes and 
                filter by AP occurences

                logical() is required to transform binary vector to a valid filter.
                %}
        indexes = 1:length(timeOfAPs);
        filter = logical(timeOfAPs(:) >= onset & timeOfAPs(:) < next_onset);
        indexesOfAPs = indexes(filter);

        iteration.APs = zeros(dataTicks,1);
        %normalize the APs, to the start of window
        APsInWindow = timeOfAPs(indexesOfAPs)-onset + 1; %the index is 1 based
        APsInWindow(APsInWindow>dataTicks)=[]; %remove aps out of ticks range
        iteration.APs(APsInWindow)=1;

        %create fixed data set
        res = NaN(dataTicks, 3);
        res(:,1) = iteration.time + onset;
        res(:,2) = iteration.stimuli;
        res(:,3) = iteration.APs;

        iteration.all = res(res(:,2) > 0 | res(:,3) > 0,:); %keep records with data
        
        %store in neuronData output
        nextIterationIndex = lastIterationIndex + length(iteration.all) + 1;
        neuronData.all(lastIterationIndex+1:nextIterationIndex-1,:) = iteration.all;
        lastIterationIndex = nextIterationIndex-1;
    end %iIteration

    %remove NaN rows
    Simulation.Neuron{iNeuron}.Data = neuronData.all(any(~isnan(neuronData.all),2),:);
end %iNeuron

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save('PreProcessed_NonRep.mat', 'Simulation');
end
