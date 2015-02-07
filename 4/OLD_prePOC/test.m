clc;
close all;

%% create the fit
for iNeuron=1:NEURONS
    curNeuron = Simulation.Neuron{iNeuron};
    data = curNeuron.CurveFitData;    


    % based on http://stackoverflow.com/questions/22792020/matlab-accumarray-weighted-mean
    %data = [10 1;30 1;20 2;30 1;20 4;20 6]
    %data = funcXAnyYData;
    %data = [0.1,10; 0.15 20; -0.5 5;0.77 2; 0.8 2;-0.44 12];

    %{

    data =

        10     1
        30     1
        20     2
        30     1
        20     4
        20     6    
    %}

    data = sortrows(data, 1);
    XVals = data(:,1); %after linear filter
    YVals = data(:,2); %raw stims
    %[funcXData, ~, Groups] = unique(data(:,1),'stable');

    BIN_SIZE=0.1;
    %BIN_SIZE=0.05; 
    %BIN_SIZE=0.001; 

    firstBin = floor(min(XVals)*10)/10;
    lastBin = ceil(max(XVals)*10)/10;
    %firstBin = min(XVals);
    %lastBin = max(XVals);
    bins = firstBin:BIN_SIZE:lastBin;
    %no zero intesity, somehow it does not equal zero, but a really small number
    %bins(abs(bins(:)*10)<BIN_SIZE)=[]; 
    [bincounts,binIndex] = histc(XVals,bins);
    %bins(bincounts(:)==0)=[]; 
    %bins = bins(unique(binIndex(:,1),'stable'));
    
    %sumByBins = accumarray(binIndex,accAPsAtLeastOne, [length(bins) 1]);
    %{
    ID =

        10
        20
        30


    Groups =

         1
         2
         2
         2
         3
         3
    %}
    %fnMean = @(ii) mean(YVals(ii));

    funcXData = bins';
    %funcXData(bincounts(:,1)==0)=[];
    %funcYMeans = accumarray(binIndex, YVals, [], @mean);
    funcYMeans = accumarray(binIndex, YVals, [length(bins) 1], @mean);
    
    %remove empty/'zero' intesity buckets
    m = [funcXData funcYMeans];
    m = m(m(:,2)~=0,:);
    funcXData = m(:,1);
    funcYMeans = m(:,2);

    %TODO: rename variables

    %figure;
    %plot(funcXData,funcYMeans,'.');

    [fitresult, gof] = createFit(funcXData,funcYMeans,iNeuron,MODE);
    curNeuron.Generator = fitresult;
    
    
    %[fitresult2, gof2] = createFitGideon(funcXData,funcYMeans);
    %[fitresult3, gof3] = createNonLinearFit(funcXData,funcYMeans);
    
    Simulation.Neuron{iNeuron} = curNeuron;
end

