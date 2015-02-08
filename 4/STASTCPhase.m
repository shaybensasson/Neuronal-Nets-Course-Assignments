close all;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)))
    clearvars -except MODE;
    load(['PreProcessed_' MODE '.mat'])
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

SAFETY_WINDOW_TO_THE_PAST_IN_STIMS = STIMS_IN_STA_WINDOW*3;
SAFETY_WINDOW_TO_THE_PAST_IN_TICKS = SAFETY_WINDOW_TO_THE_PAST_IN_STIMS * ...
            Simulation.STIMULUS_EACH_TICKS;
    
%TODO: for iNeuron=1:NEURONS
for iNeuron=2:2
    
    data = Simulation.Neuron{iNeuron}.Data;
    %get rid of data that has less history than STA safety window
    idxFirstAP = find(data(:,3)==1, 100, 'first');
    idxFirstAP(idxFirstAP<=SAFETY_WINDOW_TO_THE_PAST_IN_STIMS) = [];
    idxFirstAP = idxFirstAP(1);
    
    idxLastAP = find(data(:,3)==1, 100, 'last');
    to = idxLastAP(end)-SAFETY_WINDOW_TO_THE_PAST_IN_STIMS;
    idxLastAP(idxLastAP>to) = [];
    idxLastAP = idxLastAP(end);
    
    %mean over NaNs
    rawStimuliMean = mean(data(~isnan(data(:,2)),2));
        
    %adding index column
    data = [data (1:length(data))'];
    
    %get only AP times (from the first AP that has safety window to the
    %past)
    apsTimes = data(idxFirstAP:idxLastAP,:);
    apsTimes = apsTimes(apsTimes(:,3)==1, [1 4]);
    
    NUM_OF_APS = length(apsTimes);
    
    %accumulate all Spike Triggered stimuli for every spike
    accSTAStims = NaN(NUM_OF_APS, STIMS_IN_STA_WINDOW);
    accSTAWindow = zeros(1,STA_WINDOW_IN_TICKS+SAFETY_WINDOW_TO_THE_PAST_IN_TICKS);
    accSTAWindowCounters = zeros(1,STA_WINDOW_IN_TICKS+SAFETY_WINDOW_TO_THE_PAST_IN_TICKS);
    
    %{
	NOTE: checking must be done only for the non-rep case,
			the rep window checking and validation is not mandatory
	%}
	idxStimuliWindow = 1;
    
    for iAP=1:NUM_OF_APS %for every AP
        
        timeOfAP = apsTimes(iAP,1);
        timeOfStimuliWindow = StimTime(idxStimuliWindow);
        timeOfStimuliNextWindow = StimTime(idxStimuliWindow+1);
        
        if (timeOfAP>=timeOfStimuliNextWindow) %ap after stimuli window ends
            %go to next window
            idxStimuliWindow = idxStimuliWindow+1;
            
            timeOfStimuliWindow = StimTime(idxStimuliWindow);
            timeOfStimuliNextWindow = StimTime(idxStimuliWindow+1);
        end
                
        %get Spike Triggered stims
        idxOfAP = apsTimes(iAP,2);
        apWindowStims = data(idxOfAP-SAFETY_WINDOW_TO_THE_PAST_IN_STIMS:idxOfAP-1, [1 2]);
        apWindowStims(isnan(apWindowStims(:,2)), :)=[];
        apWindowStims = apWindowStims(find(apWindowStims(:,2), STIMS_IN_STA_WINDOW, 'last'),:);
        
        if (apWindowStims(1,1)<timeOfStimuliWindow) %discard stims before stimuli window starts
            continue;
        end
        
        accSTAStims(iAP,:) = apWindowStims(:,2)';
        
        apWindowStims(:,1) = timeOfAP - apWindowStims(:,1);
        accSTAWindow(apWindowStims(:,1)) = accSTAWindow(apWindowStims(:,1)) + apWindowStims(:,2)';
        accSTAWindowCounters(apWindowStims(:,1)) = ...
            accSTAWindowCounters(apWindowStims(:,1)) + logical(apWindowStims(:,1))';
                
    end
    
    %clear nans
    accSTAStims(all(isnan(accSTAStims),2),:)=[];
    
    accSTAWindow = accSTAWindow(1:STA_WINDOW_IN_TICKS);
    accSTAWindowCounters = accSTAWindowCounters(1:STA_WINDOW_IN_TICKS);
    
    accSTAWindow(accSTAWindowCounters==0) = NaN;
    accSTAWindowCounters(accSTAWindowCounters==0) = NaN;
    accSTAWindow = accSTAWindow./accSTAWindowCounters;
    
    %% STA creation
    %binnify the STA
    %NOTE: we can choose even smaller binSize (ex. 100 is still OK)
    BIN_SIZE = Simulation.STIMULUS_EACH_TICKS;
    bins = 1:BIN_SIZE:STA_WINDOW_IN_TICKS;
    [bincounts,binIndex] = histc(1:STA_WINDOW_IN_TICKS,bins);

    %insert any unbinned data to last bin
    unbinnedIndex = max(binIndex);
    %last bin was not created
    if (~bincounts(end))
        unbinnedIndex = unbinnedIndex+1;
    end
    binIndex(binIndex==0)=unbinnedIndex;

    %group by times and mean stim vals
    meanIgnoreNaNs = @(vector) mean(vector(~isnan(vector(:))));
    grp_mean = accumarray(binIndex', accSTAWindow', [length(bins) 1], meanIgnoreNaNs, NaN);

    STA = grp_mean';
    
    %Normalize the raw stims by the overall mean
    STA = STA - rawStimuliMean;
    
    %{
    normalize: 
        this way we unify all the filters from different neurons,
        so we could compare them.
    %}
    STA=(STA-mean(STA))./max(abs(STA));
    
    %% STC calculation
    
    %the total APs =
    % (NUM_OF_APS - the ones we discarded because they're on window edges)
    totalAPs = length(accSTAStims);
    
    %see Schwartz et al.
    stackSTA = repmat(STA,totalAPs,1);
    STC = ((accSTAStims-stackSTA)' * (accSTAStims-stackSTA)) ./(totalAPs-1);
    
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
for iNeuron=2:2
%TODO: for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    STA = Simulation.Neuron{iNeuron}.STA;
    
    STA = fliplr(STA);
    
    x = linspace(-STA_WINDOW_IN_MS, 0, length(STA));
    plot(x, STA);

    title(sprintf('STA for Neuron #%d', iNeuron));
    xlabel('Time (ms)');
    ylabel('Light levels');
 end

%% plot eVals
figure(2);

for iNeuron=2:2
%TODO: for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    evals = Simulation.Neuron{iNeuron}.EigenValues;
    
    plot(evals, 'o'); % examine eigenvalues
    
    title(sprintf('Eigenvalues (of STC) for Neuron #%d', iNeuron));
    xlabel('EigenValue/EigenVector index');
    ylabel('Variance');
 end    
    
if (SAVE_MAT_FILE)
    fprintf('Saving simulation output ...\n');
    save(['AfterSTA_' MODE '.mat'], 'Simulation');
end

%load gong 
%sound(y,Fs)
