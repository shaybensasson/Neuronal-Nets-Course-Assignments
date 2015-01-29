close all;

%works only for Rep
MODE = 'Rep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)))
    clearvars -except MODE;
    load(['PreProcessed_' MODE '.mat'])
end

NEURONS=length(Simulation.Neuron); 

STA_WINDOW_IN_MS = 500;
STA_WINDOW_IN_SEC = STA_WINDOW_IN_MS/1000;
STA_WINDOW_IN_TICKS = STA_WINDOW_IN_SEC*Simulation.TICKS_IN_SECOND;

%calculate how many stims in the STA window
STIMS_IN_STA_WINDOW = floor(STA_WINDOW_IN_TICKS/Simulation.STIMULUS_EACH_TICKS);

figure;

for iNeuron=1:NEURONS
    subplot(2,2,iNeuron);
    SAFETY_WINDOW_TO_THE_PAST = STIMS_IN_STA_WINDOW*3;
    
    data = Simulation.Neuron{iNeuron};
    %get rid of data that has less history than STA safety window
    idxFirstAP = find(data(:,3)==1, 100, 'first');
    idxFirstAP(idxFirstAP<=SAFETY_WINDOW_TO_THE_PAST) = [];
    idxFirstAP = idxFirstAP(1);
    
    idxLastAP = find(data(:,3)==1, 100, 'last');
    to = idxLastAP(end)-SAFETY_WINDOW_TO_THE_PAST;
    idxLastAP(idxLastAP>to) = [];
    idxLastAP = idxLastAP(end);
        
    %adding index column
    data = [data (1:length(data))'];
    
    %get only AP times (from the first AP that has safety window to the
    %past)
    apsTimes = data(idxFirstAP:idxLastAP,:);
    apsTimes = apsTimes(apsTimes(:,3)==1, [1 4]);
    
    NUM_OF_APS = length(apsTimes);
    
    %accumulate all Spike Triggered stimuli for every spike
    accAPs = NaN(NUM_OF_APS,1);
    accSTAStims = NaN(NUM_OF_APS, STIMS_IN_STA_WINDOW);
    totalAPs = 0;
    
    for iAP=1:NUM_OF_APS %for every AP
        
        %get Spike Triggered stims
        idxOfAP = apsTimes(iAP,2);
        windowStims = data(idxOfAP-SAFETY_WINDOW_TO_THE_PAST:idxOfAP-1, 2);
        windowStims(isnan(windowStims(:)))=[];
        windowStims = windowStims(find(windowStims, STIMS_IN_STA_WINDOW, 'last'));
        
        accSTAStims(iAP,:) = windowStims';
        
        %count the next bin (AP time + STA_WINDOW_IN_TICKS) APs
        apBin = data(idxOfAP:idxOfAP+SAFETY_WINDOW_TO_THE_PAST,[1 3]);
        apBin(apBin(:,1)<apBin(1,1)+STA_WINDOW_IN_TICKS,:)=[];
        accAPs(iAP) = length(apBin(apBin==1));
        totalAPs = totalAPs+accAPs(iAP);
    end
    
    %normalize
    accSTAStims = accSTAStims-mean(mean(accSTAStims,1));
    
    %see http://en.wikipedia.org/wiki/Spike-triggered_average
    STA = 1/totalAPs * accSTAStims'  * accAPs;

    % normalize spike-triggered average
    normSTA=(STA-mean(STA))./mean(STA);

    %plot STA
    x = linspace(-STA_WINDOW_IN_MS, 0, length(normSTA));
    plot(x, normSTA);

    title(sprintf('STA for Neuron #%d', iNeuron));
    xlabel('Time (ms)');
    ylabel('Light levels');
end

load gong 
sound(y,Fs)
        