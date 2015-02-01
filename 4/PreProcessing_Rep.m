clear all;
if (~exist('TTRep','var'))
    load('FixedData.mat')
    
    %fix diffs to be 200s instead of 100s, each time we load raw data
    StimTimeRep = StimTimeRep(1:2:end); 
end

ConstantsHeader();

SAVE_MAT_FILE = 1;

%secs in stimulus rep window
SECONDS_IN_WINDOW = 200;
TICKS_IN_SECOND = 10000;
TICKS_IN_WINDOW = TICKS_IN_SECOND * SECONDS_IN_WINDOW;
STIMULI_PER_SECOND = 1/30;
%how many ticks there are between stimuli.
STIMULUS_EACH_TICKS = TICKS_IN_SECOND*STIMULI_PER_SECOND; %333.333
STIMULI_PER_WINDOW = SECONDS_IN_WINDOW / STIMULI_PER_SECOND; %6000 stims for 200 secs

NEURONS = length(TTRep);
ITERATIONS = length(StimTimeRep); %NOTE: each ITERATION is a repetition

%window stim values by time, warming it up, to improve prefs
globalStimValues = StimulusRep(1:STIMULI_PER_WINDOW);

clear Simulation;
Simulation.Mode = CONSTANTS.MODES.REP;
Simulation.Phase = CONSTANTS.PHASES.PREPROCESSING;
Simulation.Neuron = cell(1,length(TTRep));
Simulation.StimTimeRep = StimTimeRep;
Simulation.ITERATIONS = ITERATIONS;
Simulation.SECONDS_IN_WINDOW = SECONDS_IN_WINDOW;
Simulation.TICKS_IN_WINDOW = TICKS_IN_WINDOW;
Simulation.TICKS_IN_SECOND = TICKS_IN_SECOND;
Simulation.STIMULUS_EACH_TICKS = STIMULUS_EACH_TICKS;

for iNeuron = 1:length(TTRep)
    fprintf('[N:#%i] ...\n', iNeuron);
    
    %allocate an array to store prepocessed data for neuron
    neuronData.all = NaN(STIMULI_PER_WINDOW*ITERATIONS, 3);
    lastIterationIndex = 0;

    for iIteration = 2:ITERATIONS
        %fprintf('[N:#%i] processing iteration #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 

        %actual ticks in window limited by TICKS_IN_WINDOW
        dataTicks = TICKS_IN_WINDOW;

        %Normalize stim time to start from 1        
        StimTimeRepNorm=StimTimeRep-StimTimeRep(iStart)+1;
        iteration.time= (StimTimeRepNorm(iStart):StimTimeRepNorm(iEnd)-1)'; %first 200 secs
        dataTicks = min(dataTicks,length(iteration.time)); 
        iteration.time = iteration.time(1:dataTicks); %trim

        %copy global var
        stimValues = globalStimValues(:);

        %% get stim values of the current rep window        
        iteration.stimuli = NaN(dataTicks,1);
	
        %create stim times for every stim        
        stimTimes = 1:floor(STIMULUS_EACH_TICKS):dataTicks;
        stimTimes = stimTimes(1:STIMULI_PER_WINDOW);

        iteration.stimuli(stimTimes)=stimValues;

        %% get aps of the current rep window        
        timeOfAPs=TTRep(iNeuron).sp; %time of Aps

        onset = StimTimeRep(iStart);
        next_onset = StimTimeRep(iEnd);
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
    save('PreProcessed_Rep.mat', 'Simulation');
end