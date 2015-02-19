%load('MatFiles\AfterSTA_NonRep.mat') 
MODE='NonRep';
Simulation = Sim_NonRep;
NEURONS=length(Simulation.Neuron); 
STA_WINDOW_IN_MS = Simulation.STA_WINDOW_IN_MS;

minSTAValue = inf; maxSTAValue = -inf;
minSTCFilterValue = inf; maxSTCFilterValue = -inf;

for iNeuron=1:NEURONS
    
    STA = Simulation.Neuron{iNeuron}.STA;
    
    minSTAValue = min(minSTAValue, min(STA));
    maxSTAValue = max(maxSTAValue, max(STA));
    
    Filter = Simulation.Neuron{iNeuron}.STCFilter;
    
    minSTCFilterValue = min(minSTCFilterValue, min(Filter));
    maxSTCFilterValue = max(maxSTCFilterValue, max(Filter));
end

%% plot STA
figure(1);

for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    STA = Simulation.Neuron{iNeuron}.STA;
        
    x = linspace(-STA_WINDOW_IN_MS, 0, length(STA));
    plot(x, STA);
    
    ylim([minSTAValue maxSTAValue]);
    

    title(sprintf('Neuron #%d', iNeuron));
    xlabel('Time (ms)');
    ylabel('Light levels');
    
    CreateTitleForSubplots('\bf STA');
end
 
saveas(1, ['STA_' MODE], 'png');

%% plot eVals
h = figure(2);

for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    eVals = Simulation.Neuron{iNeuron}.EigenValues;
    
    plot(eVals, 'o'); % examine eigenvalues
    
    ylim([min(eVals) max(eVals)+100]);
    
    
    title(sprintf('Neuron #%d', iNeuron));
    xlabel('EigenValue/EigenVector index');
    ylabel('Variance');
    
    CreateTitleForSubplots('\bf Eigenvalues (of STC)');
end

saveas(h, ['EigenValues_' MODE], 'png');

%% plot STC filter
h = figure(3);

for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    Filter = Simulation.Neuron{iNeuron}.STCFilter;
            
    x = linspace(-STA_WINDOW_IN_MS, 0, length(Filter));
    plot(x, Filter, 'r');
    ylim([minSTCFilterValue maxSTCFilterValue]);
   
    title(sprintf('Neuron #%d', iNeuron));
    xlabel('Time (ms)');
    ylabel('Light levels');
    
    CreateTitleForSubplots('\bf STC Filter');
end

saveas(h, ['STCFilter_' MODE], 'png');

%load gong 
%sound(y,Fs)
