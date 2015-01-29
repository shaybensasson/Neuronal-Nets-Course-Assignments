close all;

%choose Rep or NonRep
MODE = 'Rep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)))
    clearvars -except MODE;
    load(['PreProcessed_' MODE '.mat'])
end

NEURONS=length(Simulation.Neuron); 

STA_WINDOW_IN_MS = 1000;
STA_WINDOW_IN_SEC = STA_WINDOW_IN_MS/1000;
STA_WINDOW_IN_TICKS = STA_WINDOW_IN_SEC*Simulation.TICKS_IN_SECOND;

%calculate how many stims in the STA window
STIMS_IN_STA_WINDOW = floor(STA_WINDOW_IN_TICKS/Simulation.STIMULUS_EACH_TICKS);

figure;

for iNeuron=1:NEURONS
    subplot(2,2,iNeuron);
    
    data = Simulation.Neuron{iNeuron};
    %get read of data that has less history than STA window
    data = data(data(:,1)-STA_WINDOW_IN_TICKS>data(1,1), :);

    %get only AP times
    apsTimes = data(data(:,3)==1, 1);

    %get stims and their times
    stimAndTimes = data(~isnan(data(:,2)), 1:2);



    %accumulate all Spike Triggered stimuli for every spike
    accSTAData = NaN(STIMS_IN_STA_WINDOW*length(apsTimes),2);
    idxSTAChunk=0;
    NUM_OF_APS = length(apsTimes);
    REPORT_EVERY_APS = floor(NUM_OF_APS/5);
    for i=1:NUM_OF_APS %for every AP

        %get only times&stims in window
        windowStart = apsTimes(i)-STA_WINDOW_IN_TICKS;
        windowEnd = apsTimes(i);
        filter = logical(stimAndTimes(:,1)>=windowStart & stimAndTimes(:,1)<windowEnd);
        windowStims = stimAndTimes(filter, :);

        %normalize stim times
        delta = zeros(length(windowStims), 2);
        delta(:,1)=windowStart-1;
        windowStims = windowStims-delta;

        if (~mod(i, REPORT_EVERY_APS))
            fprintf('[N:#%d]apsTime #%d/%d ...\n', iNeuron, i, NUM_OF_APS);
        end

        accSTAData(idxSTAChunk+1:idxSTAChunk+length(windowStims), :) = windowStims;
        idxSTAChunk = idxSTAChunk + length(windowStims);
    end

    %group stim vals by time, and perform mean on each time stim vals
    accSTAData = sortrows(accSTAData, 1);
    Vals = accSTAData(:,2);
    [ID, ~, Groups] = unique(accSTAData(:,1),'stable');
    fnMean = @(ii) mean(Vals(ii));
    means = accumarray(Groups, 1:numel(Groups), [], fnMean);


    %average stim values across spikes, in bins of stimulus duration
    %BIN_SIZE = 10000/30;
    bins = 0:Simulation.STIMULUS_EACH_TICKS:STA_WINDOW_IN_TICKS;
    [bincounts,binIndex] = histc(1:length(means),bins);
    %binIndex(binIndex == 0) = []; %throw un-binned stims
    STA = accumarray(binIndex',means(1:length(binIndex)),[], @mean);

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
        