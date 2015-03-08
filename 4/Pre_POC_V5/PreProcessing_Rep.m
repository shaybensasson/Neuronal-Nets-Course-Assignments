clearvars -except Sim_NonRep;
if (~exist('TTRep','var'))
    load('FixedData.mat')
    
    %fix diffs to be 200s instead of 100s, each time we load raw data
    StimTimeRep = StimTimeRep(1:2:end); 
end

ConstantsHeader();

SAVE_MAT_FILE = 1;

%secs in stimulus rep trail
SECONDS_IN_TRAIL = 200;
TICKS_IN_SECOND = 10000;
TICKS_IN_TRAIL = TICKS_IN_SECOND * SECONDS_IN_TRAIL;
STIMULI_PER_SECOND = 1/30;
%how many ticks there are between stimuli.
STIMULUS_EACH_TICKS = TICKS_IN_SECOND*STIMULI_PER_SECOND; %333.333
STIMULI_PER_TRAIL = SECONDS_IN_TRAIL / STIMULI_PER_SECOND; %6000 stims for 200 secs

%throw 30 secs on the begining of each iteration, malformed data
TICKS_TO_THROW = 30 * TICKS_IN_SECOND;

NEURONS = length(TTRep);
ITERATIONS = length(StimTimeRep); %NOTE: each ITERATION is a repetition
ITERATIONS = ITERATIONS-1; %just for simplicity of loops later

clear Simulation;
Simulation.Mode = CONSTANTS.MODES.REP;
Simulation.Phase = CONSTANTS.PHASES.PREPROCESSING;
Simulation.Neuron = cell(1,length(TTRep));
Simulation.StimTimeRep = StimTimeRep;
Simulation.StimulusRep = StimulusRep;
Simulation.TT = TTRep;
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
   
    for iIteration = 1:ITERATIONS
        %fprintf('[N:#%i] iter #%i/#%i ...\n', iNeuron, iIteration, ITERATIONS);

        iStart = iIteration; 
        iEnd = iIteration+1;

        %actual ticks in trail limited by TICKS_IN_TRAIL
        dataTicks = TICKS_IN_TRAIL;

        iteration.time= (StimTimeNonRep(iStart):StimTimeNonRep(iEnd)-1)'; %first 200 secs
        dataTicks = min(dataTicks,length(iteration.time)); 

        %% get stim values of the current rep trail        
        stimValues = StimulusRep(1:STIMULI_PER_TRAIL,:);

        %Smoothen stim values on dataTicks
        stimValues = repmat(stimValues,1,ceil(STIMULUS_EACH_TICKS));
        stimValues = stimValues';
        stimValues = stimValues(:);
        
        stimValues(ceil(STIMULUS_EACH_TICKS)*2:ceil(STIMULUS_EACH_TICKS)*3:end)=[];
        stimValues(ceil(STIMULUS_EACH_TICKS):(ceil(STIMULUS_EACH_TICKS)-1)+ceil(STIMULUS_EACH_TICKS)*2:end)=[];
        
        %because of flooring we might have less data
        dataTicks = min(dataTicks,length(stimValues)); 
        
        %trim data collected so far
        iteration.time = iteration.time(1:dataTicks);
        iteration.stimuli=stimValues(1:dataTicks);

        %% get aps of the current rep trail        
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

        iteration.APs = NaN(dataTicks,1);

        %convert times to indexes
        indexesOfAPsInTrail = timeOfAPs(indexesOfAPs)-onset + 1; %the index is 1 based
        indexesOfAPsInTrail(indexesOfAPsInTrail>dataTicks)=[]; %remove aps out of ticks range
        iteration.APs(indexesOfAPsInTrail)=1;

        %create fixed data set
        iteration.all = NaN(dataTicks, 3);
        iteration.all(:,1) = iteration.time;
        iteration.all(:,2) = iteration.stimuli;
        iteration.all(:,3) = iteration.APs;

        %throw initial ticks
        iteration.all(1:TICKS_TO_THROW,:) = [];
        
        Simulation.Neuron{iNeuron}.Iteration{iStart} = iteration.all;
    end %iIteration
end %iNeuron

Sim_Rep = Simulation;

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save('MatFiles\PreProcessed_Rep.mat', 'Sim_Rep', '-v7.3');
end

load gong 
sound(y,Fs)