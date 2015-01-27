clc;
close all;
clear all;

if ((~exist('StimTime','var')) || (~exist('TT','var')))
    load('SpikeTimeFFFONOFF.mat')
end

INTERACTIVE = 1;

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

figureEx('Custom', 'Maximize');
%figure

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

%% plot neurons data
for iNeuron=1:NUM_OF_NEURONS
    
    subplot(NUM_OF_NEURONS+1,1,iNeuron+1);
    title(sprintf('Neuron #%d', iNeuron));
    
    hold on
    
    ylabel('# of Iteration')
    xlabel('Time from stimulus onset (s)')
    axis([0 STIMULUS_WINDOW 0 ITERATIONS])
    set(gca,'XTickLabel',sprintf('%1.1f|',0:0.2:2))
    
    timeOfAPs=TT(iNeuron).sp;
    indexesPerNeuron = 1:length(timeOfAPs);
    
    for iIteration=1:ITERATIONS
        currentOff = StimTime(iIteration);
        currentOn = currentOff+DELAY_TO_ON;
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
            %normalize the APs (to the start of the current window)
            normalizedAPs = timeOfAPs(indexesOfAPs)-currentOff;
            
            indexes = 1:length(normalizedAPs);
            
            indexesOn = indexes(normalizedAPs(:) >= currentOn & ...
                    normalizedAPs(:) < nextOff);
            scatter(normalizedAPs(indexesOn), ...
                ones(length(normalizedAPs(indexesOn)),1)*iIteration,...
                1,'k');
            
            indexesOff = indexes(normalizedAPs(:) < currentOn);
            scatter(normalizedAPs(indexesOff), ...
                ones(length(normalizedAPs(indexesOff)),1)*iIteration, ...
                1,'k');
        end
                
        if (INTERACTIVE)
            if (~mod(iIteration,10))
                drawnow;
            end
        else
            if (~mod(iIteration,10))
                fprintf('[n=%d,i=%d/%d] plotting ...\n', ...
                    iNeuron, iIteration, ITERATIONS);
            end
        end
        
        
    end %for ITERATIONS
    
    
    maxy = ITERATIONS;
    axis([0 STIMULUS_WINDOW 0 maxy]);
    med = round(maxy/2);
    set(gca,'YTick',[0 med maxy]);
    
end %for iNeuron

        
beep('on');

