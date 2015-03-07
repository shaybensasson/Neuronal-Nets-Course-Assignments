%{
res = repmat(STA', 1, ceil(Simulation.STIMULUS_EACH_TICKS));
res = res'
res = res(:)';
%}

%{
%take STA_WINDOW_IN_TICKS before AP
windowOnSet = timeOfAP-1 -STA_WINDOW_IN_TICKS;
%apWindowStims(:,1) = apWindowStims(:,1)-windowOnSet;

%because we took some safety stims, let's take the last relevant ones
apWindowStims = apWindowStims(find(apWindowStims(:,1),Simulation.STIMS_IN_STA_WINDOW+1,'last'), :);

%let's smoothen the STA from Simulation.STIMS_IN_STA_WINDOW+1 into STA_WINDOW_IN_TICKS
stims = apWindowStims(:,2);
stims = repmat(stims, 1, floor(Simulation.STIMULUS_EACH_TICKS));
stims = stims';
stims = stims(:)';
start = length(stims)-STA_WINDOW_IN_TICKS;
stims = stims(start:start+STA_WINDOW_IN_TICKS-1);
%}
%{
stimsOnSet = timeOfAP-STA_WINDOW_IN_TICKS-ceil(Simulation.STIMULUS_EACH_TICKS)

stims = apWindowStims(apWindowStims(:,1)>=stimsOnSet, 2);
res = repmat(stims, 1, ceil(Simulation.STIMULUS_EACH_TICKS));
res = res';
res = res(:)';
start = length(res)-STA_WINDOW_IN_TICKS;
res = res(start:start+STA_WINDOW_IN_TICKS-1);
%}

%stimValues = repmat(stimValues,1,floor(STIMULUS_EACH_TICKS));
stimValues = repmat(stimValues,1,2);
stimValues = stimValues';
stimValues = stimValues(:);