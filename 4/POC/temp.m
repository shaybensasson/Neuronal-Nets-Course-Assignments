NUM_OF_BIN_SIZES = numel(Simulation.PSTH_BIN_SIZES);
NUM_OF_NEURONS = length(Simulation.Neuron);
result = NaN(NUM_OF_BIN_SIZES*NUM_OF_NEURONS, 5);
idx = 1;
for iBinSize = 1:numel(Simulation.PSTH_BIN_SIZES)
    for iNeuron = 1:length(Simulation.Neuron)
        data = Simulation.Neuron{iNeuron}.PSTH{iBinSize};
        result(idx, :) = [ceil(data.BIN_SIZE) iNeuron ...
            data.SSEk data.SSEg ...
            data.SSEk-data.SSEg];
        idx = idx+1;
    end
end