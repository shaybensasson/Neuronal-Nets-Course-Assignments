close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

load(['MatFiles\PreProcessed_' MODE '.mat'])
    
switch MODE
    case 'Rep'
        Simulation = Sim_Rep;
        clearvars Sim_Rep;
        StimTime = Simulation.StimTimeRep;
        
    case 'NonRep'
        Simulation = Sim_NonRep;
        clearvars Sim_NonRep;
        StimTime = Simulation.StimTimeNonRep;
        
    otherwise
        ME = MException('noSuchMODE', ...
            'no such MODE is found!');
        throw(ME)
end

SAVE_MAT_FILE = 1;
CREATE_STC = 1;

NEURONS=length(Simulation.Neuron); 

STA_WINDOW_IN_MS = 1000;
STA_WINDOW_IN_SEC = STA_WINDOW_IN_MS/1000;
STA_WINDOW_IN_TICKS = STA_WINDOW_IN_SEC*Simulation.TICKS_IN_SECOND;

Simulation.Phase = CONSTANTS.PHASES.STASTC;
Simulation.STA_WINDOW_IN_MS = STA_WINDOW_IN_MS;
Simulation.STA_WINDOW_IN_SEC = STA_WINDOW_IN_SEC;
Simulation.STA_WINDOW_IN_TICKS = STA_WINDOW_IN_TICKS;

minSTAValue = inf; maxSTAValue = -inf;
minSTCFilterValue = inf; maxSTCFilterValue = -inf;
        
for iNeuron=1:NEURONS
   
    fprintf('[N:#%i] ...\n', iNeuron);
    
    %accumulate all Spike Triggered stimuli for every spike
    STAStims = NaN(length(Simulation.TT(iNeuron).sp), STA_WINDOW_IN_TICKS);
    
    STAPerTrail = zeros(Simulation.ITERATIONS,STA_WINDOW_IN_TICKS);
    countTrailAPs = 0;
    STA = zeros(1,STA_WINDOW_IN_TICKS);
    
    %the total APs =
    % (NUM_OF_APS - the ones we discarded because they're on trail edges)
    totalAPsCounted = 0;
           
    for iIteration=1:Simulation.ITERATIONS
        
        data = Simulation.Neuron{iNeuron}.Iteration{iIteration};
        
        %start from AP that has enough history
        idxFirstAP = find(data(STA_WINDOW_IN_TICKS+1:end,3)==1, 1, 'first') + STA_WINDOW_IN_TICKS;
        
        %adding index column
        %discard warning we have 4 neurons
        data = [data (1:length(data))'];

        %get only AP times 
        apsTimes = data(idxFirstAP:end,:);
        apsTimes = apsTimes(apsTimes(:,3)==1, [1 4]);

        NUM_OF_APS = length(apsTimes);
        
        for iAP=1:NUM_OF_APS %for every AP

            timeOfAP = apsTimes(iAP,1);

            %get Spike Triggered stims
            idxOfAP = apsTimes(iAP,2);

            apStims = data(idxOfAP-STA_WINDOW_IN_TICKS:idxOfAP-1, [1 2]);

            stims = apStims(:,2)';

            totalAPsCounted = totalAPsCounted+1;
            STAStims(totalAPsCounted,:) = stims;

            STAPerTrail(iIteration,:) =  ...
                STAPerTrail(iIteration,:) + stims;
            countTrailAPs = countTrailAPs+1;
            
        end %iAP
        
        STA = STA + STAPerTrail(iIteration,:)./countTrailAPs;
    
    end %iIteration
    
    if (CREATE_STC) %pref improvement
        %get the data physically used
        STAStims = STAStims(1:totalAPsCounted, :);
    end
    
    %{
    normalize: 
        this way we unify all the filters from different neurons,
        so we could compare them.
    
        NOTE: checked without normalize, the linear fit looks worse
    %}
    STA=(STA-mean(STA))./max(abs(STA));
    
    %use for ploting
    minSTAValue = min(minSTAValue, min(STA));
    maxSTAValue = max(maxSTAValue, max(STA));
    
    
    Simulation.Neuron{iNeuron}.STA = STA;

    %% STC calculation
    if (CREATE_STC)

        fprintf('STC calculation ...\n');
        
        %FUTURE compare vs pillow: 
        %see Schwartz et al.
        %stackSTA = repmat(STA,totalAPs,1);
        %STC = ((accSTAStims-stackSTA)' * (accSTAStims-stackSTA)) ./(totalAPs-1);

        %see
        %http://pillowlab.cps.utexas.edu/teaching/CompNeuro10/slides/slides07_STC_LNPmodel.pdf,
        %slide 43
        STC = STAStims' * STAStims ./(totalAPsCounted);

        [V,D] = eig(STC);
        %{ 
        produces matrices of eigenvalues (D, diag(D) are eigenvalues) 
        and eigenvectors (V, its columns are the eigenvectors of A) of matrix A, 
        so that A*V = V*D.
        %}

        eVals = diag(D); %just as SVD, these are variances
        eVects = V; 

        %look for the most variance Ev
        [~,idx] = max(eVals);

        Filter = eVects(:, idx)';
        Filter = normalize(Filter, 1, 1);

        %use for ploting
        minSTCFilterValue = min(minSTCFilterValue, min(Filter));
        maxSTCFilterValue = max(maxSTCFilterValue, max(Filter));

        Simulation.Neuron{iNeuron}.EigenValues = eVals;
        Simulation.Neuron{iNeuron}.STCFilter = Filter;
    end
    
end %iNeuron

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

if (CREATE_STC)
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
end

%% save
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
    save(['MatFiles\AfterSTA_' MODE '.mat'], ['Sim_' MODE], '-v7.3');
end

load gong 
sound(y,Fs)
