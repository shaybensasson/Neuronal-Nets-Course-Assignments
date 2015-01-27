clc;
close all;
clear all;

if ((~exist('StimTime','var')) || (~exist('TT','var')))
    load('SpikeTimeFFFONOFF.mat')
end

%{
TT(i): contains data recorded from neuron i.
TT(i).sp: the times of the APs
a. times are in 1/10^4 second
b. from the onset of the experiment

StimTime: the time that the light was set to OFF, after the experiment
        onset.
StimTime(i) = the time that the light was set to OFF on the i^th iteration
StimTime(i)+10^4 = time that the light was set back to ON on the i^th iteration

%}

NUM_OF_NEURONS=length(TT);
DELAY_TO_ON=10000;
STIMULUS_WINDOW = 20000;
ITERATIONS = length(StimTime); 

BIN_SIZE=100;
NUM_OF_BINS = STIMULUS_WINDOW/BIN_SIZE;
bins = 0:BIN_SIZE:STIMULUS_WINDOW;
indexesOfStimulusTime = 1:STIMULUS_WINDOW;

figureEx('Custom', 'Maximize');
%figure;

%% plot Light stimulus bars
subplot(NUM_OF_NEURONS+1,1,1);
onOffData = zeros(STIMULUS_WINDOW,1);
onOffData(1:DELAY_TO_ON)=1;
h1=bar(onOffData,'k');
hold on;
onOffData(1:DELAY_TO_ON)=0;
onOffData(DELAY_TO_ON+1:STIMULUS_WINDOW)=1;
h2=bar(onOffData);
set(h2,'FaceColor','y','EdgeColor','y');
axis([0 STIMULUS_WINDOW 0 1])
set(gca,'XTickLabel',sprintf('%1.1f|',0:0.2:2));
set(gca,'YTickLabel','');
title('Light stimulus on/off');

maxHz = 0; %Max Hz for graphs

%% plot neurons data
for iNeuron=1:NUM_OF_NEURONS
    
    subplot(NUM_OF_NEURONS+1,1,iNeuron+1);
    title(sprintf('Neuron #%d (Bin size: %d*10^{-4})', iNeuron, BIN_SIZE));
    
    hold on
    ylabel('Firing rate (Hz)')
    xlabel('Time from stimulus onset (s)')
    set(gca,'XTickLabel',sprintf('%1.1f|',0:0.2:2))
    
    timeOfAPs=TT(iNeuron).sp;
    indexesPerNeuron = 1:length(timeOfAPs);
    accAPs = zeros(STIMULUS_WINDOW, 1);
        
    for iIteration=1:ITERATIONS
        currentOff = StimTime(iIteration);
        nextOn_onset = currentOff+DELAY_TO_ON;
        nextOff = currentOff+STIMULUS_WINDOW;
        
        %{
        Change from time scale into discrete indexes and 
        filter by AP occurences
        
        logical() is required to transform binary vector to a valid filter.
        %}
        filter = logical(timeOfAPs(:,1) >= currentOff & timeOfAPs(:,1) < nextOff);
        indexesOfAPs = indexesPerNeuron(filter);
                
        %we have any APs
        if(~isempty(indexesOfAPs))
            %normalize the APs
            normalizedAPs = timeOfAPs(indexesOfAPs)-currentOff;
            
            curAPs = zeros(STIMULUS_WINDOW, 1);
            curAPs(normalizedAPs) = 1;
            accAPs = accAPs + curAPs;
        end
    end %for ITERATIONS
    
    %plot historam for neuron
    filter = accAPs(:,1) > 0;
    spikeTimes = indexesOfStimulusTime(filter)';
    accAPsAtLeastOne = accAPs(filter);
    
    %{
    fprintf('[n=%d] sum accAPs=%d ...\n', ...
        iNeuron, sum(accAPs));
    %}
    
    [bincounts,binIndex] = histc(spikeTimes,bins);
    sumByBins = accumarray(binIndex,accAPsAtLeastOne, [length(bins) 1]);
    bar(bins,sumByBins,'histc');
    maxHz = max(max(sumByBins), maxHz);
    axis([0 STIMULUS_WINDOW 0 maxHz]);
    med = round(maxHz/2);
    set(gca,'YTick',[0 med maxHz]);
    set(gca,'YTickLabel',['  ';num2str(med);num2str(maxHz)]);

    
end

        
beep('on');

