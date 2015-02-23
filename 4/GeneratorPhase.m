close all;
clearvars -except Sim_Rep Sim_NonRep;
ConstantsHeader();

%choose Rep or NonRep
MODE = 'NonRep';

load(['MatFiles\AfterLinearFilter_' MODE '.mat'])
    
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

NEURONS = length(Simulation.Neuron);
             
%store for later usage
Simulation.Phase = CONSTANTS.PHASES.GENERATOR;

for iBinSize=1:numel(Simulation.RATE_BIN_SIZES)
    curBinSize = Simulation.RATE_BIN_SIZES(iBinSize);
    
    title = sprintf('Generator (%s), Bin size: %.2f, Using %s Filter', ...
        MODE, curBinSize, iif(Simulation.UsingSTA,'STA','STC'));
    
    hf = figure(iBinSize);
    hf.Name = title;
    
    for iNeuron=1:NEURONS
        subplot(2,2,iNeuron);
        curNeuron = Simulation.Neuron{iNeuron};
        fprintf('[Bsz,N:#%d,#%d] ...\n', iBinSize, iNeuron);

        %get binned psth and stims after K
        rateData = curNeuron.Rate{iBinSize}.Data;
        psth = rateData(:,2);
        stimsAfterLinearFilter = rateData(:,3);
        
        %% create the fit
        curveFitData = [stimsAfterLinearFilter psth];

        curveFitData = sortrows(curveFitData, 1);
        XVals = curveFitData(:,1); %after linear filter
        YVals = curveFitData(:,2); %psth

        %NOTE: kept for diag later
        %plot(XVals, YVals, 'o'); %we get a shiftet gausian
        
        FIT_BIN_SIZE=0.1;

        %create nice looking numbers
        firstBin = floor(min(XVals)*10)/10;
        lastBin = ceil(max(XVals)*10)/10;

        %put into bins
        bins = firstBin:FIT_BIN_SIZE:lastBin;
        [bincounts,binIndex] = histc(XVals,bins);
        
        %insert any unbinned data to last bin (we're gonna remove it later)
        %binIndex(binIndex==0)=max(binIndex);

        %group by mean
        funcXData = bins';
        funcYMeans = accumarray(binIndex, YVals, [length(bins) 1], @mean, NaN);

        %remove empty/'zero' intesity buckets
        m = [funcXData funcYMeans bincounts];
        %m = m(2:end-1, :); %throw first and last bins, really few stims there
        m = m(~isnan(m(:,2)),:);
        
        funcXData = m(:,1);
        funcYMeans = m(:,2);
        bincounts = m(:,3);
                
        %create the 'function' using Interpolant Linear Fit
        [fitresult, gof] = createFit(funcXData,funcYMeans,bincounts,iNeuron,MODE);
        curNeuron.Rate{iBinSize}.Generator = fitresult;
                
        
        %% apply the generator
        NORMALIZE = 1; NORMALIZE_BY_MAX=1;
        stimsAfterGenerator = normalize( ...
            fitresult(stimsAfterLinearFilter), ...
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
    
    CreateTitleForSubplots(['\bf ' title]);
    
    r = 150; %pixels pre inch
    set(hf, 'PaperUnits', 'inches');
    set(hf, 'PaperPosition', [0 0 2880 1620]/r); %x_width=10cm y_width=15cm
    
    saveas(hf, ['Generator_' MODE '_BinSize_' ...
        sprintf('%d', floor(curBinSize)) ...
        '_UsingSTA_' num2str(Simulation.UsingSTA)], 'png');
    
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
