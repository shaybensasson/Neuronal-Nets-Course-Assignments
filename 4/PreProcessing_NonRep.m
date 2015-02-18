clearvars -except Sim_Rep;
if (~exist('TTNonRep','var'))
    load('FixedData.mat')
end

ConstantsHeader();

SAVE_MAT_FILE = 1;

%secs in stimulus non-rep trail
SECONDS_IN_TRAIL = 100;
TICKS_IN_SECOND = 10000;
TICKS_IN_TRAIL = TICKS_IN_SECOND * SECONDS_IN_TRAIL;
STIMULI_PER_SECOND = 1/30;
%how many ticks there are between stimuli.
STIMULUS_EACH_TICKS = TICKS_IN_SECOND*STIMULI_PER_SECOND; %333.333
STIMULI_PER_TRAIL = SECONDS_IN_TRAIL / STIMULI_PER_SECOND; %3000 stims for 100 secs

%throw 30 secs on the begining of each iteration, malformed data
TICKS_TO_THROW = 30 * TICKS_IN_SECOND;

NEURONS = length(TTNonRep);
ITERATIONS = length(StimTimeNonRep); %NOTE: each ITERATION is different bulk (non-repeating trail)
ITERATIONS = ITERATIONS-2; %we ignore the last 2 chunks because it has partial stimuli values

clear Simulation;
Simulation.Mode = CONSTANTS.MODES.NONREP;
Simulation.Phase = CONSTANTS.PHASES.PREPROCESSING;
Simulation.Neuron = cell(1,length(TTNonRep));
Simulation.StimTimeNonRep = StimTimeNonRep;
Simulation.ITERATIONS = ITERATIONS;
Simulation.SECONDS_IN_TRAIL = SECONDS_IN_TRAIL;
Simulation.TICKS_IN_TRAIL = TICKS_IN_TRAIL;
Simulation.TICKS_IN_SECOND = TICKS_IN_SECOND;
Simulation.STIMULI_PER_SECOND = STIMULI_PER_SECOND;
Simulation.STIMULUS_EACH_TICKS = STIMULUS_EACH_TICKS;
Simulation.STIMULI_PER_TRAIL = STIMULI_PER_TRAIL;
Simulation.TICKS_TO_THROW = TICKS_TO_THROW;

        
for iNeuron = 1:length(TTNonRep)
    
    fprintf('[N:#%i] ...\n', iNeuron);
    
    lastIterationSimuliIndex = 0;
    
    for iIteration = 2:ITERATIONS
        %fprintf('[N:#%i] iter #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

        iEnd = iIteration;
        iStart = iEnd-1;
        
        %actual ticks in trail limited by TICKS_IN_TRAIL
        dataTicks = TICKS_IN_TRAIL;
        
        %Normalize stim time to start from 1
        StimTimeNonRepNorm=StimTimeNonRep-StimTimeNonRep(iStart)+1;
        iteration.time= (StimTimeNonRepNorm(iStart):StimTimeNonRepNorm(iEnd)-1)'; %first 100 secs
        dataTicks = min(dataTicks,length(iteration.time)); 
        iteration.time = iteration.time(1:dataTicks); %trim
        
        %% get stim values of the current non-rep trail        
        nextIterationSimuliIndex = lastIterationSimuliIndex + STIMULI_PER_TRAIL + 1;
        stimValues = StimulusNonRep(lastIterationSimuliIndex+1:nextIterationSimuliIndex-1,:);
        lastIterationSimuliIndex = nextIterationSimuliIndex-1;
        
        %Smoothen stim values on dataTicks
        %TODO: we might try with floor + 1
        stimValues = repmat(stimValues,1,ceil(STIMULUS_EACH_TICKS));
        stimValues = stimValues';
        stimValues = stimValues(:);
        iteration.stimuli=stimValues(1:dataTicks);

        %% get aps of the current non-rep trail        
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

        iteration.APs = NaN(dataTicks,1);
        %normalize the APs, to the start of trail
        indexesOfAPsInTrail = timeOfAPs(indexesOfAPs)-onset + 1; %the index is 1 based
        indexesOfAPsInTrail(indexesOfAPsInTrail>dataTicks)=[]; %remove aps out of ticks range
        iteration.APs(indexesOfAPsInTrail)=1;

        %create fixed data set
        iteration.all = NaN(dataTicks, 3);
        iteration.all(:,1) = iteration.time + onset -1;
        iteration.all(:,2) = iteration.stimuli;
        iteration.all(:,3) = iteration.APs;
        
        %throw initial ticks
        iteration.all(1:TICKS_TO_THROW,:) = [];
        
        Simulation.Neuron{iNeuron}.Iteration{iStart} = iteration.all;
    end %iIteration
end %iNeuron

Sim_NonRep = Simulation;

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save('MatFiles\PreProcessed_NonRep.mat', 'Sim_NonRep');
end
