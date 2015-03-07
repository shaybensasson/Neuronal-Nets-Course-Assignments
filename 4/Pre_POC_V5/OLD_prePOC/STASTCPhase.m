close all;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'Rep';

if (~exist('Simulation','var') || (~strcmp(Simulation.Mode,MODE)) || ...
        Simulation.Phase < CONSTANTS.PHASES.PSTH)
    clearvars -except MODE;
    load(['AfterPSTH_' MODE '.mat'])
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

for iNeuron=1:NEURONS
    SAFETY_WINDOW_TO_THE_PAST = STIMS_IN_STA_WINDOW*3;
    
    data = Simulation.Neuron{iNeuron}.Data;
    %get rid of data that has less history than STA safety window
    idxFirstAP = find(data(:,3)==1, 100, 'first');
    idxFirstAP(idxFirstAP<=SAFETY_WINDOW_TO_THE_PAST) = [];
    idxFirstAP = idxFirstAP(1);
    
    idxLastAP = find(data(:,3)==1, 100, 'last');
    to = idxLastAP(end)-SAFETY_WINDOW_TO_THE_PAST;
    idxLastAP(idxLastAP>to) = [];
    idxLastAP = idxLastAP(end);
        
    %adding index column
    data = [data (1:length(data))'];
    
    %get only AP times (from the first AP that has safety window to the
    %past)
    apsTimes = data(idxFirstAP:idxLastAP,:);
    apsTimes = apsTimes(apsTimes(:,3)==1, [1 4]);
    
    NUM_OF_APS = length(apsTimes);
    
    %accumulate all Spike Triggered stimuli for every spike
    accAPs = NaN(NUM_OF_APS,1);
    accSTAStims = NaN(NUM_OF_APS, STIMS_IN_STA_WINDOW);
    totalAPs = 0;
    
    %{
	NOTE: checking must be done only for the non-rep case,
			the rep window checking and validation is not mandatory
	%}
	idxStimuliWindow = 1;
    
    for iAP=1:NUM_OF_APS %for every AP
        
        timeOfAP = apsTimes(iAP,1);
        timeOfStimuliWindow = StimTime(idxStimuliWindow);
        timeOfStimuliNextWindow = StimTime(idxStimuliWindow+1);
                
        %get Spike Triggered stims
        idxOfAP = apsTimes(iAP,2);
        apWindowStims = data(idxOfAP-SAFETY_WINDOW_TO_THE_PAST:idxOfAP-1, [1 2]);
        apWindowStims(isnan(apWindowStims(:,2)), :)=[];
        apWindowStims = apWindowStims(find(apWindowStims(:,2), STIMS_IN_STA_WINDOW, 'last'),:);
        
        if (apWindowStims(1,1)<timeOfStimuliWindow) %discard stims before stimuli window starts
            continue;
        end
        
        
        %count the next bin (AP time + STA_WINDOW_IN_TICKS) APs
        apBin = data(idxOfAP:idxOfAP+SAFETY_WINDOW_TO_THE_PAST,[1 3]);
        apBin(apBin(:,1)<apBin(1,1)+STA_WINDOW_IN_TICKS,:)=[];
        
        if (apBin(end,1)>=timeOfStimuliNextWindow) %discard stims after stimuli window ends
            %go to next window if we're after it
            idxStimuliWindow = idxStimuliWindow+1;
            continue;
        end
        
        accSTAStims(iAP,:) = apWindowStims(:,2)';
        accAPs(iAP) = length(apBin(apBin==1));
        totalAPs = totalAPs+accAPs(iAP);
    end
    
	%clear nans
    accSTAStims(all(isnan(accSTAStims),2),:)=[];
    accAPs(all(isnan(accAPs),2),:)=[];
    
    %{
    normalize: 
        this way we unify all the filters from different neurons,
        so we could compare them.
    %}
    accSTAStims = accSTAStims-mean(mean(accSTAStims,1));
    
    %see http://en.wikipedia.org/wiki/Spike-triggered_average
    STA = 1/totalAPs * accSTAStims'  * accAPs;

    % normalize spike-triggered average
    STA=(STA-mean(STA))./max(abs(STA));
    
    
    STC = accSTAStims'*(accSTAStims.*repmat(accAPs,1,STIMS_IN_STA_WINDOW)) ...
        / (totalAPs-1) - STA*STA'*totalAPs/(totalAPs-1);

    Simulation.Neuron{iNeuron}.STA = STA;
    Simulation.Neuron{iNeuron}.STC = STC;
    
    [V,D] = eig(STC);
    %{ 
    produces matrices of eigenvalues (D, diag(D) are eigenvalues) 
    and eigenvectors (V, its columns are the eigenvectors of A) of matrix A, 
    so that A*V = V*D.
    %}
    
    evals = abs(diag(D)); %just as SVD, these are variances
    Simulation.Neuron{iNeuron}.EigenValues = evals;
    evects = V; 
    Simulation.Neuron{iNeuron}.EigenVectors = evects;
end %iNeuron

%% plot STA
figure(1);
for iNeuron=1:NEURONS
    subplot(2,2, iNeuron);
    STA = Simulation.Neuron{iNeuron}.STA;
    
    x = linspace(-STA_WINDOW_IN_MS, 0, length(STA));
    plot(x, STA);

    title(sprintf('STA for Neuron #%d', iNeuron));
    xlabel('Time (ms)');
    ylabel('Light levels');
 end

%% plot eVals
figure(2);
for iNeuron=1:NEURONS
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
