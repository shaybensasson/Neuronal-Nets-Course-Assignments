if (~exist('StimTimeRep','var'))
    %load('fixedData.mat');
    load('Data.mat');
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

figureEx('units','normalized','outerposition',[0 0 1 1])

%normalized into 200 diffs, currently 100 diffs
StimTimeRepEvery200 = StimTimeRep(1:2:end);
%StimTimeRepEvery200(2)-StimTimeRepEvery200(1) ~= 200*10^4

%sampling rate is 1/10000, we multiply by 3 so it'd be 1/30000
time = 3;

