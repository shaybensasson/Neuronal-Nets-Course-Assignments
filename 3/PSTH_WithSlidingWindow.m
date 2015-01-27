clc;
close all;
clear all;

if ((~exist('StimTime','var')) || (~exist('TT','var')))
    load('SpikeTimeFFFONOFF.mat')
end

%{
TT(i): contains data recorded from neuron i.
TT(i).sp: the times of the APs
a. times are in 1/10^5 second
b. from the onset of the experiment

StimTime: the time that the light was set to OFF, after the experiment
        onset.
StimTime(i) = the time that the light was set to OFF on the i^th iteration
StimTime(i)+10^5 = time that the light was set back to ON on the i^th iteration

%}

NUM_OF_NEURONS=length(TT);
DELAY_TO_ON=10000;
STIMULUS_WINDOW = 20000;
ITERATIONS = length(StimTime); 

GAUSSIAN_WINDOW = 0; %Gaussian or sliding window (MA)
GAUSSIAN_SIGMA = 1;
SLIDING_WINDOW = 5; %MA assumes w is odd, gaus shall be pair

%plots a gaussian with params
if (GAUSSIAN_WINDOW)
    x = -round(SLIDING_WINDOW/2):.1:round(SLIDING_WINDOW/2);
    norm = normpdf(x,0,GAUSSIAN_SIGMA);
    figure;
    plot(x,norm)
end
    

BIN_SIZE=100;

NUM_OF_BINS = STIMULUS_WINDOW/BIN_SIZE;
bins = 0:BIN_SIZE:STIMULUS_WINDOW;
indexesOfStimulusTime = 1:STIMULUS_WINDOW;

figureEx('Custom', 'Maximize');

%% plot Light stimulus bars
subplot(NUM_OF_NEURONS+1,1,1);
onOffData = zeros(STIMULUS_WINDOW,1);
onOffData(1:DELAY_TO_ON)=1;
h1=bar(onOffData,'y');
hold on;
onOffData(1:DELAY_TO_ON)=0;
onOffData(DELAY_TO_ON+1:STIMULUS_WINDOW)=1;
h2=bar(onOffData);
set(h2,'FaceColor','k','EdgeColor','k');
axis([0 STIMULUS_WINDOW 0 1])
set(gca,'XTickLabel',sprintf('%1.1f|',0:0.2:2));
set(gca,'YTickLabel','');
title('Light stimulus on/off');

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
        On_onset = currentOff-DELAY_TO_ON;
        nextOn_onset = currentOff+DELAY_TO_ON;
        
        %{
        Change from time scale into discrete indexes and 
        filter by AP occurences
        
        logical() is required to transform binary vector to a valid filter.
        %}
        filter = logical(timeOfAPs(:,1) >= On_onset & timeOfAPs(:,1) < nextOn_onset);
        indexesOfAPs = indexesPerNeuron(filter);
        
        %normalize: as timespan started from the current On onset
        normalizedOff_onset = currentOff-On_onset;
                
        
        %we have any APs
        if(~isempty(indexesOfAPs))
            %normalize the APs
            normalizedAPs = timeOfAPs(indexesOfAPs)-On_onset;
            
            curAPs = zeros(STIMULUS_WINDOW, 1);
            curAPs(normalizedAPs) = 1;
            accAPs = accAPs + curAPs;
        end
                
        %move to next chunk
        On_onset = nextOn_onset;
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
    
    if (GAUSSIAN_WINDOW) 
        sumByBinsFiltered = slidingGaussian(sumByBins, SLIDING_WINDOW, ...
            GAUSSIAN_SIGMA);
    else %MovingAverage
        sumByBinsFiltered = movingAverage(sumByBins, SLIDING_WINDOW);
    end
    h1 = bar(bins,sumByBins,'histc');
    set(h1,'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0.60 0.60 0.60], ...
        'FaceAlpha',0.50);
    hold on;
    h2 = bar(bins,sumByBinsFiltered,'histc');
    set(h2,'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','FaceAlpha',1);

    maxy = round(max(sumByBins));
    axis([0 STIMULUS_WINDOW 0 maxy]);
    med = round(maxy/2);
    set(gca,'YTick',[0 med maxy]);

    
end

        
beep('on');

