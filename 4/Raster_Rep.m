if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,'Rep')))
    clear all;
    load('PreProcessed_Rep.mat')
end

INTERACTIVE = 0;

NEURONS = length(Simulation.Neuron);
ITERATIONS = Simulation.ITERATIONS;
StimTimeRep = Simulation.StimTimeRep;
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;

maxWindowLength = 0;

%% plot neurons data
for iNeuron=1:NEURONS
    subplot(NEURONS,1,iNeuron);
    title(sprintf('Neuron #%d', iNeuron));
    
    hold on
    
    ylabel('# of Iteration')
    xlabel('Time from stimulus onset (s)')
    axis([0 TICKS_IN_WINDOW 0 ITERATIONS])
    set(gca,'XTickLabel',sprintf('%1.1f|',0:20:SECONDS_IN_WINDOW));
    
    for iIteration=2:ITERATIONS
        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 

        times = Simulation.Neuron{iNeuron}(:,1);
        
        onset = StimTimeRep(iStart);
        next_onset = StimTimeRep(iEnd);
        filter = logical(times(:) >= onset & times(:) < next_onset);
        windowData = Simulation.Neuron{iNeuron}(filter, :);
        
        delta = zeros(length(windowData), 3);
        delta(:,1) = onset;
        windowData=windowData-delta; %normalize
        
        actualWindowLength = StimTimeRep(iEnd)-StimTimeRep(iStart)+1;
        maxWindowLength = max(maxWindowLength, actualWindowLength);
        
        scatter(windowData(:,1), ...
                windowData(:,3)*iIteration,...
                1,'k');

                
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
    axis([0 TICKS_IN_WINDOW 0 maxy]);
    med = round(maxy/2);
    set(gca,'YTick',[0 med maxy]);
    
end %for iNeuron

        
beep('on');

