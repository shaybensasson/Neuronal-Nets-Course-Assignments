close all;
clear all;
clc;

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

figure;
%figureEx('Custom', 'Maximize');

%sampling rate is 1/10000, we multiply by 3 so it'd be 1/30000, because
%10000/30 is not an integer
TIME_FACTOR = 3; 
SAMPLING_TIME_RESOLUTION = 10000;
TICKS_IN_SECOND=TIME_FACTOR*SAMPLING_TIME_RESOLUTION; %ex. 30000 per sec
REPETITION_WINDOW_IN_SECS = 200; %in secs

NUM_OF_NEURONS=4; 
REPETITION_WINDOW=REPETITION_WINDOW_IN_SECS*TICKS_IN_SECOND; %200*3*10^4=6*10^6

STIMULUS_INTERVAL_IN_SECS = 1/30; %stimulus appears every 1/30 s
STIMULUS_INTERVAL = TICKS_IN_SECOND*STIMULUS_INTERVAL_IN_SECS; % in ticks

%normalized into 200 diffs, currently 100 diffs
StimTimeRepFixed = StimTimeRep(1:2:end);
StimTimeRepFixed = StimTimeRepFixed*TIME_FACTOR;
%StimTimeRepFixed(2)-StimTimeRepFixed(1) ~= 200*10^4  *TIME_FACTOR

NUM_OF_STIMULI_VALUES = 6000;

% how long do we look back at the stimulus from the ap?
STA_WINDOW_IN_MSECS = 1000;  %500(expected:15) 400(expected:12)
MSEC_IN_SEC = 1000; %msecs in 1 sec
STA_WINDOW=(TICKS_IN_SECOND/MSEC_IN_SEC)*STA_WINDOW_IN_MSECS; %12000 equals 400 m.s.
EXPECTED_STIMULI_COUNT = STA_WINDOW/STIMULUS_INTERVAL; %ex. 12

%normalize
%StimulusRepNormalized=StimulusRep-mean(StimulusRep);
StimulusRepNormalized=StimulusRep;

stimuliTimes = zeros(NUM_OF_STIMULI_VALUES,1);
for iNeuron=1:NUM_OF_NEURONS
    subplot(2,2,iNeuron);
    
    staSum=zeros(EXPECTED_STIMULI_COUNT,1);
        
    timeOfAPs=TTRep(iNeuron).sp*TIME_FACTOR;
    %prune every ap that doesn't have full window of STA before it
    timeOfAPs(timeOfAPs < StimTimeRepFixed(1)+STA_WINDOW) = []; 
    
    TOTAL_APS = length(timeOfAPs);
    
    for iAPTime=1:TOTAL_APS
        
        
        curAPTime = timeOfAPs(iAPTime);
        
        %we look back STA_WINDOW time
        staBegin = curAPTime-STA_WINDOW;
        
        %lets find the last stimulus that occured before the AP
        iClosestStimTime = find(StimTimeRepFixed(:)<curAPTime,1,'last');
        if(staBegin<StimTimeRepFixed(iClosestStimTime))
            %the closest stimulus we found is inside the window, 
            %we shall discard the ap (too fast for our STA window).
            %We'd like that STA will begin inside the window, 
            %so we can find intersection stims and the ap
            %{
            fprintf( ...
                        ['staBegin<=StimTimeRepFixed(iClosestStimTime)' ...
                        ' N(%d),APTime(%d)\n'], ...
                           iNeuron, iAPTime);
            %}
        else %STA begins inside stim time repetition window
            
            %define a window backwards in time STA_WINDOW long
            %we take -1, because of causality - stim causes ap
            currentWindow = curAPTime-STA_WINDOW:curAPTime-1;
            
            %get the stim repetition window starting from the closest stim time 
            % (jumping in STIMULUS_INTERVAL, occurs every 1/30 sec)
            stimuliTimes = ... 
                (StimTimeRepFixed(iClosestStimTime):STIMULUS_INTERVAL: ...
                (StimTimeRepFixed(iClosestStimTime)+REPETITION_WINDOW)-1)';
            
            indexesOfStimTimesAtSTA = find(ismember(stimuliTimes,currentWindow));
                        
            countSimuli = length(indexesOfStimTimesAtSTA);
            if (countSimuli<EXPECTED_STIMULI_COUNT)
                    fprintf( ...
                        ['NOTE: stimuliAtSTA != %d: %d.' ...
                        ' N(%d),APTime(%d)\n'], ...
                           EXPECTED_STIMULI_COUNT, ...
                           countSimuli, iNeuron, iAPTime);
            end

            if(countSimuli==(EXPECTED_STIMULI_COUNT-1)) 
                %we expect EXPECTED_STIMULI_COUNT, but we substract 1
                % at the end of 'currentWindow', due to causality.
                %ex. so instead of 12000/1000=12 will have 11 stimuli
                % to solve this we add an earlier datapoint (index)
                temp=zeros(EXPECTED_STIMULI_COUNT,1);
                temp(2:EXPECTED_STIMULI_COUNT)=indexesOfStimTimesAtSTA;
                temp(1)=temp(2)-1;
                indexesOfStimTimesAtSTA=temp;
            end

            %get the stim values of the STA
            stimuliValuesAtSTA=StimulusRepNormalized(indexesOfStimTimesAtSTA);
            
            %{
            if(length(stimuliValuesAtSTA)~=EXPECTED_STIMULI_COUNT)
                fprintf( ...
                    ['WARNING: length(stimuliValuesAtSTA)<%d>~=EXPECTED_STIMULI_COUNT<%d>.' ...
                    ' N(%d),APTime(%d)\n'], ...
                    length(stimuliValuesAtSTA), EXPECTED_STIMULI_COUNT, iNeuron, iAPTime);
            end
            %}
                
            %accumulate all values in STA window
            staSum = staSum+stimuliValuesAtSTA;
        
        end
         
    end
    
    
    %{
    fprintf( ...
                    ['Discarded %d Aps.' ...
                        ' N(%d)\n'], ...
                        REPETITION_WINDOW/STIMULUS_INTERVAL, iNeuron);
    %}
    
    staSum = staSum ./ TOTAL_APS;
        
    x = linspace(-STA_WINDOW_IN_MSECS, 0, EXPECTED_STIMULI_COUNT);
    plot(x, staSum , '-', 'LineWidth',1.5);
    hold on;
    %scatter(x, staSum, 40, 'r');
    title(sprintf('STA for Neuron #%d', iNeuron));
    xlabel('Time (ms)');
    ylabel('Light levels');
end