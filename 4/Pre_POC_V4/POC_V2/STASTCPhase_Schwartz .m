close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

load(['FiringRate_' MODE '.mat'])
    
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

NEURONS=length(Simulation.Neuron); 

STA_WINDOW_IN_MS = 1000;
STA_WINDOW_IN_SEC = STA_WINDOW_IN_MS/1000;
STA_WINDOW_IN_TICKS = STA_WINDOW_IN_SEC*Simulation.TICKS_IN_SECOND;

%calculate how many stims in the STA window
STIMS_IN_STA_WINDOW = floor(STA_WINDOW_IN_TICKS/Simulation.STIMULUS_EACH_TICKS);

Simulation.Phase = CONSTANTS.PHASES.STASTC;
Simulation.STA_WINDOW_IN_MS = STA_WINDOW_IN_MS;
Simulation.STA_WINDOW_IN_SEC = STA_WINDOW_IN_TICKS;
Simulation.STA_WINDOW_IN_TICKS = STA_WINDOW_IN_TICKS;
Simulation.STIMS_IN_STA_WINDOW = STIMS_IN_STA_WINDOW;

Simulation.PSTH_BIN_SIZES = [Simulation.STIMULUS_EACH_TICKS; ... %sampling freq
             Simulation.STIMULUS_EACH_TICKS*3; ... 100 msec          
             Simulation.STIMULUS_EACH_TICKS*5; ... ~150 msec
             Simulation.STIMULUS_EACH_TICKS*10; ... ~330 msec
             Simulation.STIMULUS_EACH_TICKS*15; ... ~500 msec
             Simulation.STIMULUS_EACH_TICKS*30]; %1000 msec = 1 sec 

SAFETY_WINDOW_TO_THE_PAST_IN_STIMS = STIMS_IN_STA_WINDOW*3;
SAFETY_WINDOW_TO_THE_PAST_IN_TICKS = SAFETY_WINDOW_TO_THE_PAST_IN_STIMS * ...
            Simulation.STIMULUS_EACH_TICKS;
    
%for iNeuron=2:2
for iNeuron=1:NEURONS
   
    fprintf('[N:#%i] ...\n', iNeuron);
    
    data = Simulation.Neuron{iNeuron}.RawData;
    
    %mean over NaNs
    rawStimuliMean = mean(data(~isnan(data(:,2)),2));
        
    %adding index column
    %discard warning we have 4 neurons
    data = [data (1:length(data))'];
    
    %get only stim data
    dataForStimsOnly = data(~isnan(data(:,2)),[1 2 4]);
    
    totalBins = dataForStimsOnly(end,1)-dataForStimsOnly(1,1);
    totalBins = ceil(totalBins/Simulation.STIMULUS_EACH_TICKS/STIMS_IN_STA_WINDOW);
    
    totalStims = length(dataForStimsOnly);
    
    %accumulate all Spike Triggered stimuli for each bin
    accSTAStims = NaN(totalBins, STIMS_IN_STA_WINDOW);
    totalAps = 0; %total aps counted
    totalActualBins = 0; %actual bins omitting the discarded ones
    
    %{
	NOTE: checking must be done only for the non-rep case,
			the rep window checking and validation is not mandatory
	%}
	idxStimuliWindow = 1;
    
    for iBin=2:totalBins %for every bin
        
        timeOfStimuliWindow = StimTime(idxStimuliWindow);
        timeOfStimuliNextWindow = StimTime(idxStimuliWindow+1);
        
        idxFirstStim = ((iBin-2)*STIMS_IN_STA_WINDOW)+1;
        idxLastStim = idxFirstStim+STIMS_IN_STA_WINDOW-1;
        
        %get stims on bin
        idxFirstStimOnNextBin = idxLastStim+1;
        idxLastStimOnNextBin = idxFirstStimOnNextBin+STIMS_IN_STA_WINDOW-1;
        
        if (idxLastStimOnNextBin > totalStims)
            %we've gone thru all stimuli
            break;
        end
            
        stimsOnBin = dataForStimsOnly(idxFirstStim:idxLastStim, 2)'; 
        
        %The APs that the stims in the prev bin triggered
        apWindowFromOnData = dataForStimsOnly(idxLastStim,3)+1;
        apWindowToOnData = dataForStimsOnly(idxLastStimOnNextBin,3);
        
        
        timeOfLastAP = data(apWindowToOnData, 1);
        
        if (timeOfLastAP>timeOfStimuliNextWindow) %ap after stimuli window ends
            %go to next window
            idxStimuliWindow = idxStimuliWindow+1;
            
            %discard this bin
            continue;
        end
        
        totalActualBins = totalActualBins+1;
        
        totalApsOnNextBin = nansum(data(apWindowFromOnData:apWindowToOnData, 3));
        accSTAStims(totalActualBins,:) = stimsOnBin * totalApsOnNextBin;
        totalAps = totalAps+totalApsOnNextBin;
        
    end
    
    %clear nans
    accSTAStims = accSTAStims(1:totalActualBins, :);
    
    STA = sum(accSTAStims)*(1/totalBins);
    
    %Normalize by the overall mean
    STA = STA - rawStimuliMean;
    
    %{
    normalize: 
        this way we unify all the filters from different neurons,
        so we could compare them.
    %}
    STA=(STA-mean(STA))./max(abs(STA));
    
    %% STC calculation
    
    
    %see Schwartz et al.
    stackSTA = repmat(STA,length(accSTAStims),1);
    STC = ((accSTAStims-stackSTA)' * (accSTAStims-stackSTA)) ./(totalAps-1);
    
    Simulation.Neuron{iNeuron}.STA = STA;
    Simulation.Neuron{iNeuron}.STC = STC;
    
    [V,D] = eig(STC);
    %{ 
    produces matrices of eigenvalues (D, diag(D) are eigenvalues) 
    and eigenvectors (V, its columns are the eigenvectors of A) of matrix A, 
    so that A*V = V*D.
    %}
    
    evals = diag(D); %just as SVD, these are variances
    Simulation.Neuron{iNeuron}.EigenValues = evals;
    evects = V; 
    Simulation.Neuron{iNeuron}.EigenVectors = evects;
end %iNeuron

%% plot STA
figure(1);

%for iNeuron=2:2
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    STA = Simulation.Neuron{iNeuron}.STA;
    
    %NOTE: we flip before plotting
    STA = fliplr(STA);
    
    x = linspace(-STA_WINDOW_IN_MS, 0, length(STA));
    plot(x, STA);

    title(sprintf('STA for Neuron #%d', iNeuron));
    xlabel('Time (ms)');
    ylabel('Light levels');
 end

%% plot eVals
figure(2);

%for iNeuron=2:2
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    evals = Simulation.Neuron{iNeuron}.EigenValues;
    
    plot(evals, 'o'); % examine eigenvalues
    text(1:length(evals), evals+1,num2str(evals, '%.4f'), 'Rotation', 45);
    
    title(sprintf('Eigenvalues (of STC) for Neuron #%d', iNeuron));
    xlabel('EigenValue/EigenVector index');
    ylabel('Variance');
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

%% save
if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['AfterSTA_' MODE '.mat'], ['Sim_' MODE]);
end

%load gong 
%sound(y,Fs)
