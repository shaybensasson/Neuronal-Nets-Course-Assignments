if (~exist('TTNonRep','var'))
    load('FixedData.mat')
end

%Normalize stim time to start from 1
StimTimeNonRepNorm=StimTimeNonRep-StimTimeNonRep(1);
%ExperimentTicks=1:max(StimTimeNonRepNorm);
%ExperimentTicks(ExperimentTicks==StimTimeNonRepNorm)
%indexes = 1:length(StimTimeNonRepNorm);
%ind = indexes(logical(StimTimeNonRepNorm(:)));


SECONDS_IN_WINDOW = 100;
TICKS_IN_SECOND = 10000;
TICKS_IN_WINDOW = TICKS_IN_SECOND * SECONDS_IN_WINDOW;
STIMULI_PER_SECOND = 1/30;
STIMULUS_EACH_TICKS = floor(TICKS_IN_SECOND*STIMULI_PER_SECOND); %333
STIMULI_PER_WINDOW = SECONDS_IN_WINDOW / STIMULI_PER_SECOND; %3000 stims for 100 secs

simulation.time= (StimTimeNonRepNorm(1)+1:1:StimTimeNonRepNorm(2)+1)'; %first 100 secs

actualWindowLength = StimTimeNonRep(2)-StimTimeNonRep(1)+1;

%window stim values by time
stimValues= StimulusNonRep(1:STIMULI_PER_WINDOW);
stimValues = repmat(stimValues,2,1);

simulation.stimuli = NaN(actualWindowLength,1);
stimTimes = 1:STIMULUS_EACH_TICKS:actualWindowLength;
simulation.stimuli(stimTimes)=stimValues(1:length(stimTimes));

onset = StimTimeNonRep(1)+1;
next_onset = StimTimeNonRep(2)+1;

timeOfAPs=TTNonRep(1).sp; %time of Aps

%{
        Change from time scale into discrete indexes and 
        filter by AP occurences
        
        logical() is required to transform binary vector to a valid filter.
        %}
indexes = 1:length(timeOfAPs);
filter = logical(timeOfAPs(:,1) >= onset & timeOfAPs(:,1) < next_onset);
indexesOfAPs = indexes(filter);

%we have any APs
%if(~isempty(indexesOfAPs))
    
simulation.APs = zeros(actualWindowLength,1);
%normalize the APs, to the start of window
APsInWindow = timeOfAPs(indexesOfAPs)-onset;
simulation.APs(APsInWindow)=1;

%% finally, get
all = zeros(actualWindowLength, 3);
all(:,1) = simulation.time;
all(:,2) = simulation.stimuli;
all(:,3) = simulation.APs;

simulation.all = all(all(:,2) > 0 | all(:,3) > 0,:);


