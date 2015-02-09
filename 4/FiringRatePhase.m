close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

load(['PreProcessed_' MODE '.mat'])
    
switch MODE
    case 'Rep'
        Simulation = Sim_Rep;
        clearvars Sim_Rep;
    case 'NonRep'
        Simulation = Sim_NonRep;
        clearvars Sim_NonRep;
    otherwise
        ME = MException('noSuchMODE', ...
            'no such MODE is found!');
        throw(ME)
end

SAVE_MAT_FILE = 1;
NEURONS = length(Simulation.Neuron);

%store for later usage
Simulation.Phase = CONSTANTS.PHASES.GENERATOR;
Simulation.RATE_BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq
             Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec          
             Simulation.STIMULUS_EACH_TICKS*5; ... ~150 msec
             Simulation.STIMULUS_EACH_TICKS*10; ... ~330 msec
             Simulation.STIMULUS_EACH_TICKS*15; ... ~500 msec
             Simulation.STIMULUS_EACH_TICKS*30]; %1000 msec = 1 sec 
         
for iBinSize=1:numel(Simulation.RATE_BIN_SIZES)
    curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    
    for iNeuron = 1:NEURONS
        fprintf('[Bsz,N:#%d,#%d] ...\n', iBinSize, iNeuron);
        
        curNeuron = Simulation.Neuron{iNeuron};
        
        times = curNeuron.RawData(:,1);
        stimValues = curNeuron.RawData(:,2);
        aps = curNeuron.RawData(:,3);
        
        %group raw data times by bins
        bins = times(1):curBinSize:times(end);
        [~,binIndex] = histc(times,bins);

        %insert any unbinned data to last bin
        binIndex(binIndex==0)=max(binIndex);

        %group by times and sum aps
        grp_psth = accumarray(binIndex, aps, [length(bins) 1], @nansum, NaN);

        rateData = [bins' grp_psth];
        
        %remove NaNs
        rateData(isnan(rateData(:, 2)), :) = [];

        %throw bins before first ap
        firstApTime = curNeuron.RawData(~isnan(aps) & aps>0, 1);
        firstApTime = firstApTime(1);

        %throw bins before first ap
        %FUTURE: why do we check for psth == 0?
        rateData(rateData(:,1)<firstApTime & rateData(:,2) == 0, :) = [];

        %throw last bin, it has garbage unbinned data (it's acc is NaN)
        rateData = rateData(1:end-1,:); 
        psth = rateData(:, 2);

        NORMALIZE = 1;
        NORMALIZE_BY_MAX = 1;
        psth = normalize(psth, NORMALIZE, NORMALIZE_BY_MAX);
        
        curNeuron.Rate{iBinSize}.BinSize = curBinSize;
        curNeuron.Rate{iBinSize}.Data = [rateData(:,1) psth];
        
        Simulation.Neuron{iNeuron} = curNeuron;
    end
end


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
    save(['FiringRate_' MODE '.mat'], ['Sim_' MODE]);
end
