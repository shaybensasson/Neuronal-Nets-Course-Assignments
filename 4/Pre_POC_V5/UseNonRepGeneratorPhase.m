close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%MODE is alway Rep here
MODE = 'Rep';

load('MatFiles\AfterGenerator_NonRep.mat')

load('MatFiles\AfterLinearFilter_Rep.mat')
Simulation = Sim_Rep;
clearvars Sim_Rep;
StimTime = Simulation.StimTimeRep;

SAVE_MAT_FILE = 1;

%Wether to use STA or STC?
UsingSTA = Sim_NonRep.UsingSTA;

NEURONS = length(Simulation.Neuron);
             
%store for later usage
Simulation.Phase = CONSTANTS.PHASES.GENERATOR;

for iBinSize=1:numel(Simulation.RATE_BIN_SIZES)
    curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    
    for iNeuron=1:NEURONS
        curNeuron = Simulation.Neuron{iNeuron};
        fprintf('[Bsz,N:#%d,#%d] ...\n', iBinSize, iNeuron);

        %get binned psth and stims after K
        rateData = curNeuron.Rate{iBinSize}.Data;
        psth = rateData(:,2);
        stimsAfterLinearFilter = rateData(:,3);
        
        %use appropriate generator
        generator = Sim_NonRep.Neuron{iNeuron}.Rate{iBinSize}.Generator;
        
        %% apply the generator
        NORMALIZE = 1; NORMALIZE_BY_MAX=1;
        stimsAfterGenerator = normalize( ...
            generator(stimsAfterLinearFilter), ...
            NORMALIZE, NORMALIZE_BY_MAX);
        
        rateData(:,4) = stimsAfterGenerator;
        curNeuron.Rate{iBinSize}.Data = rateData;
        
        %% create statistics
                
        %see http://www.mathworks.com/matlabcentral/answers/104189-calculate-sum-of-square-error
        SSEk = norm(psth-stimsAfterLinearFilter,2)^2; %lower is less err
        SSEg = norm(psth-stimsAfterGenerator,2)^2; %lower is less err
                
        curNeuron.Rate{iBinSize}.SSEk = SSEk;
        curNeuron.Rate{iBinSize}.SSEg = SSEg;
        
        Simulation.Neuron{iNeuron} = curNeuron;    
    end %for iNeuron
end %for iBinSize

switch MODE
    case 'Rep'
        Sim_Rep = Simulation;
    case 'NonRep'
        Sim_NonRep = Simulation;
    otherwise
        ME = MException('noSuchMODE', ...
            'no such MODE is found!');
        throw(ME)
end


if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['MatFiles\AfterGenerator_' MODE '.mat'], ['Sim_' MODE], '-v7.3');
end
        
beep('on');
