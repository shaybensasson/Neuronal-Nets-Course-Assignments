close all;

ConstantsHeader();

%choose Rep or NonRep
MODE = 'Rep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)) || ...
        Simulation.Phase < CONSTANTS.PHASES.STASTC)
    clearvars -except MODE;
    load(['AfterSTA_' MODE '.mat'])
    ConstantsHeader();
end

switch MODE
    case 'Rep'
        StimTime = Simulation.StimTimeRep;
    case 'NonRep'
        StimTime = Simulation.StimTimeNonRep;
    otherwise
        ME = MException('STA:noSuchMODE', ...
            'no such MODE is found!');
        throw(ME)
end

SAVE_MAT_FILE = 1;

NEURONS = length(Simulation.Neuron);
ITERATIONS = Simulation.ITERATIONS;
SECONDS_IN_WINDOW = Simulation.SECONDS_IN_WINDOW;
TICKS_IN_WINDOW = Simulation.TICKS_IN_WINDOW;
TICKS_IN_SECOND = Simulation.TICKS_IN_SECOND;

STA_WINDOW_IN_TICKS = Simulation.STA_WINDOW_IN_TICKS;
STIMS_IN_STA_WINDOW = Simulation.STIMS_IN_STA_WINDOW;

%throw 0-padded convolved values, see conv doc for more info;
STIMS_TO_THROW_AFTER_CONV = ceil(STIMS_IN_STA_WINDOW/2);
%duration in seconds to display when ploting multiple rates
SECONDS_OF_RATE_TO_DISPLAY = 10;

%store for later usage
Simulation.Phase = CONSTANTS.PHASES.LINEARFILTER;
Simulation.STIMS_TO_THROW_AFTER_CONV = STIMS_TO_THROW_AFTER_CONV;
Simulation.SECONDS_OF_RATE_TO_DISPLAY = SECONDS_OF_RATE_TO_DISPLAY;
        
%% apply a linear filter of STA
for iNeuron=1:NEURONS
    curNeuron = Simulation.Neuron{iNeuron};
    STA = curNeuron.STA;
    
    %create column for filtered data
    Simulation.Neuron{iNeuron}.Data(:,4) = ...
        NaN(length(Simulation.Neuron{iNeuron}.Data),1);
        
    for iIteration=2:ITERATIONS
        %Normalize stim time to start from 1
        iEnd = iIteration;
        iStart = iEnd-1; 
 
        times = curNeuron.Data(:,1);
        stimValues = curNeuron.Data(:,2);
        
        onset = StimTime(iStart);
        next_onset = StimTime(iEnd);
        
        %filter window of stimValues
        filter = logical(times(:) >= onset & times(:) < next_onset ...
            & ~isnan(stimValues(:)));
        stimValuesOnWindow = curNeuron.Data(filter, 2);
                
        %TODO: we might normalize stims before conv
        
        %{ 
        linear filter the stimVals by the STA filter
        we flip the STA so the convolution will be correct
        
        see http://en.wikipedia.org/wiki/Cross-correlation
        %}
        stimsAfterLinearFilter = conv(stimValuesOnWindow, flipud(STA'), 'valid');

        %added NaNs to 0-padded convolved values, see conv doc for more info;
        %FUTURE: we might have problems here with different sizes of STA
        stimsAfterLinearFilter = padarray(stimsAfterLinearFilter,...
            [STIMS_TO_THROW_AFTER_CONV-1 0], NaN, 'pre');
        stimsAfterLinearFilter = padarray(stimsAfterLinearFilter,...
            [STIMS_TO_THROW_AFTER_CONV 0], NaN, 'post');
        
        curNeuron.Data(filter,4) = stimsAfterLinearFilter;
    end %for ITERATIONS
    
    %NOTE: somehow the filter inserts zeros on data that 
    %  is not included in the query, so we replace them with NaNs
    curNeuron.Data(curNeuron.Data(:,4) == 0,4) = NaN;

    %remove conv padded rows
    %filter = (~isnan(curNeuron.Data(:,2))) & (isnan(curNeuron.Data(:,4)));
    filter = (curNeuron.Data(:,3) == 0) & (isnan(curNeuron.Data(:,4)));
    curNeuron.Data(filter,:)=[]; %remove padded
        
    Simulation.Neuron{iNeuron} = curNeuron;
    
end %for iNeuron

if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['AfterLinearFilter_' MODE '.mat'], 'Simulation');
end
        
beep('on');

