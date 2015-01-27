tic;
clc;
close all;
clear all;

NodesPerRowCol=20;%number of nodes in column/row
NodesPerPattern = NodesPerRowCol^2; %number of nodes in a net/pattern
%NodesPerPattern = 20;

SYNC_UPDATE = 0; %sync or async update

INITIAL_NOISE_PERCENT=0; %precentage of noise %0
NOISE_LEVELS = 1;
NOISE_STEP = 0.1;
curNoise = INITIAL_NOISE_PERCENT;

USE_NATURAL_IMAGES = 0;
BW_CONVERTION_BINARY_TH=0.4;

%Random or determintic
RANDOM_NOISE=1;
MAX_ITERATIONS = 5000;

PLOT_DIFF = 0; %plots each inferenced pattern diff
SAVE_PLOTS = 1; %whether to save plots, if 1, plots are invisible

%third question
DIVERT=0;
%shift size
MAX_DIVERT_SIZE=1;
curDivert=MAX_DIVERT_SIZE;

totalPatterns=110; %patterns stored in net

minNumOfPatterns=1;
%minNumOfPatterns=totalPatterns;

%num of iters until correct inference across all pattern
numOfIterationsUntilCorrectInf = zeros(NOISE_LEVELS,totalPatterns);
numOfCorrectInferences = zeros(NOISE_LEVELS,totalPatterns);
numOfIterationsUntilInf = zeros(NOISE_LEVELS,totalPatterns);


ERROR_TOLERANCE = 0; %amount of error in inference we tolerate

if USE_NATURAL_IMAGES
    fprintf('loading images ...\n');
    storedPatterns = getNaturalImages(NodesPerRowCol, totalPatterns,BW_CONVERTION_BINARY_TH);
else
    storedPatterns = cell(totalPatterns,1);
    for patternNum=1:totalPatterns
        storedPatterns{patternNum} = GenerateNetwork(NodesPerRowCol);
    end
end

firstRun=1;
while ((DIVERT==1&&curDivert<=MAX_DIVERT_SIZE) || firstRun==1)
    firstRun=0;
    for curNoiseLevel=1:NOISE_LEVELS
        noisyNodesNumber = ceil(NodesPerPattern*curNoise); %number of noisy nodes in a net
        
        for curNumOfPatterns=minNumOfPatterns:totalPatterns
            if(DIVERT==1)
                fprintf('[cPs:%d,cDiv:%d] Infering patterns ...\n', ...
                        curNumOfPatterns, curDivert);
            else
                fprintf('[cPs:%d,cNoi:%d] Infering patterns ...\n', ...
                        curNumOfPatterns, curNoiseLevel);
            end
            
            permPatterns = randperm(curNumOfPatterns);
            for t=1:curNumOfPatterns
                %inputPattern=round(rand(1,1)*(curNumOfPatterns-1))+1;
                inputPattern=permPatterns(t);
                %{
                fprintf('[cPs:%d,cNoi:%d] Infering pattern #%d ...\n', ...
                    curNumOfPatterns, curNoiseLevel, inputPattern);
                %}
                                
                %TODO: we might use vector multiplication
                storedPatternsWeights= cell(curNumOfPatterns,1);
                for patternNum=1:curNumOfPatterns
                    storedPatternsWeights{patternNum} = GetWeights(storedPatterns{patternNum});
                end
                
                %Sum all weights
                sumWeights=zeros(NodesPerPattern,NodesPerPattern);
                for patternNum=1:curNumOfPatterns
                    sumWeights = sumWeights+storedPatternsWeights{patternNum};
                end
                
                %Set diagnoal to zero
                sumWeights(1:NodesPerPattern+1:NodesPerPattern*NodesPerPattern) = 0;
                sumWeights = (1/NodesPerPattern)*sumWeights;
                
                %add noise
                biasedPatterns = cell(curNumOfPatterns,1);
                
                for patternNum=1:curNumOfPatterns
                    biasedPatterns{patternNum} = storedPatterns{patternNum}(:);
                    if (curNoise~=0 && ~DIVERT)
                        if(RANDOM_NOISE==1)
                            p = randperm(NodesPerPattern);
                            indexesToFlip = p(1:noisyNodesNumber);
                            biasedPatterns{patternNum}(indexesToFlip) = -1*biasedPatterns{patternNum}(indexesToFlip);
                        else
                            biasedPatterns{patternNum}(1:noisyNodesNumber) = -1*biasedPatterns{patternNum}(1:noisyNodesNumber);
                        end
                    end
                end
                
                inferedNetwork = biasedPatterns{inputPattern}(:)';
                if(DIVERT==1)
                    inferedNetwork = reshape(inferedNetwork',NodesPerRowCol,NodesPerRowCol);
                    %move 0 vertically and curDivert horiz
                    inferedNetwork = circshift(inferedNetwork, [0, curDivert]);
                    inferedNetwork = inferedNetwork(:)';
                    
                    %here we store the diverted pattern back in the
                    %biasedPatterns, so we could see it when plotting
                    biasedPatterns{inputPattern}(:) = inferedNetwork';
                end
                
                %synchornious updating
                if (SYNC_UPDATE == 1)
                    %synchornious updating
                    %inferedNetwork = W3*net3(:);
                    inferedNetwork = (2*(sumWeights*inferedNetwork' >= 0) - 1)'; %update rule
                    numOfIterationsUntilInf(curNoiseLevel,curNumOfPatterns)= numOfIterationsUntilInf(curNoiseLevel,curNumOfPatterns)+1;
                    
                    
                    correctLevel = GetCorrectPrecentage(storedPatterns{inputPattern}(:)',inferedNetwork);
                    if(correctLevel==1 )
                        numOfIterationsUntilCorrectInf(curNoiseLevel,curNumOfPatterns)= numOfIterationsUntilCorrectInf(curNoiseLevel,curNumOfPatterns)+1;
                        numOfCorrectInferences(curNoiseLevel,curNumOfPatterns)= numOfCorrectInferences(curNoiseLevel,curNumOfPatterns)+1;
                    end
                    
                    if (PLOT_DIFF)
                        if (SAVE_PLOTS)
                            h=figure('Visible','off');
                        else
                            h=figure;
                        end
                        subplot(4,1,1)
                        imagesc(storedPatterns{inputPattern})
                        title(['allNets[' num2str(inputPattern) ']']);

                        %figure
                        subplot(4,1,2)
                        imagesc(reshape(biasedPatterns{inputPattern},NodesPerRowCol,NodesPerRowCol))
                        title(['biasedPatterns[' num2str(inputPattern) ']' ...
                            iif(DIVERT, [', Divert : ' num2str(curDivert)], ...
                                [', Noise : ' num2str(curNoise)]) ...
                        ', SYNC UPDATE']);


                        %figure
                        subplot(4,1,3)
                        inferedNetworkForDisp = reshape(inferedNetwork',NodesPerRowCol,NodesPerRowCol);
                        imagesc(inferedNetworkForDisp)
                        title(['inferedNetwork' ...
                            ', SYNC UPDATE' ...
                            ', CorrectLvl: ' num2str(correctLevel)]);

                        subplot(4,1,4)
                        imagesc(storedPatterns{inputPattern}-inferedNetworkForDisp)
                        title(['DiffFromAll[' num2str(inputPattern) ']' ...
                            ', SYNC UPDATE' ...
                            ', CorrectLvl: ' num2str(correctLevel)]);
                        
                        if (SAVE_PLOTS)
                            saveas(h, ['Sync' ...
                                '_Pat' num2str(inputPattern) ...
                                iif(DIVERT, ['_Divert' num2str(curDivert)], ...
                                    ['_Noise' num2str(curNoise)]) ...
                                '~' datestr(now,'dd_mm_yy_HH_MM_SS_FFF') ...
                                '.jpg']);
                        end
                    end %if PLOT_DIFF
                else %asynchornious updating
                    
                    maxCorrectnessReached=0;
                    curIteration=1;
                    while (curIteration<MAX_ITERATIONS) && (maxCorrectnessReached==0)
                        patternNum = round(rand(1,1)*(NodesPerPattern-1))+1;

                        %all links/ksharim that involve the i neuron - sumWeights(:,i);
                        %sumweight is nxn
                        %taking the pattern/inferedNet and does * the
                        %relation of a neuron to the rest of the neurons,
                        %across patterns
                        %the patterns are stored in the sumWeights
                        %rule from wikipedia
                        %inferedNetwork is the S
                        
                        
                        nodeOnNextIteration = inferedNetwork*sumWeights(:,patternNum);
                        
                        %we get a scalar, that if positive they are on the
                        %same direction so don't change it
                        if nodeOnNextIteration>0
                            nodeOnNextIteration=1;
                        else
                            nodeOnNextIteration=-1;
                        end
                                                
                        inferedNetwork(patternNum)=nodeOnNextIteration;
                        numOfIterationsUntilInf(curNoiseLevel,curNumOfPatterns)= numOfIterationsUntilInf(curNoiseLevel,curNumOfPatterns)+1;
                        
                        correctLevel = GetCorrectPrecentage(storedPatterns{inputPattern}(:)',inferedNetwork);
                        if(correctLevel>=1-ERROR_TOLERANCE )
                            maxCorrectnessReached=1;
                            
                            numOfIterationsUntilCorrectInf(curNoiseLevel,curNumOfPatterns)= numOfIterationsUntilCorrectInf(curNoiseLevel,curNumOfPatterns)+curIteration;
                            numOfCorrectInferences(curNoiseLevel,curNumOfPatterns)= numOfCorrectInferences(curNoiseLevel,curNumOfPatterns)+1;
                        end
                        
                        %NOTE: NOT TO BE USED WITH INTENSIVE RUNS
                        if (PLOT_DIFF)
                            if (maxCorrectnessReached == 1 || curIteration == MAX_ITERATIONS-1) %disp last

                                if (SAVE_PLOTS)
                                    h=figure('Visible','off');
                                else
                                    h=figure;
                                end
                                
                                subplot(4,1,1)
                                imagesc(storedPatterns{inputPattern})
                                title(['allNets[' num2str(inputPattern) ']']);

                                %figure
                                subplot(4,1,2)
                                imagesc(reshape(biasedPatterns{inputPattern},NodesPerRowCol,NodesPerRowCol))
                                title(['biasedPatterns[' num2str(inputPattern) ']' ...
                                    iif(DIVERT, [', Divert : ' num2str(curDivert)], ...
                                        [', Noise : ' num2str(curNoise)]) ...
                                    ', Iter #' num2str(curIteration)]);


                                %figure
                                subplot(4,1,3)
                                inferedNetworkForDisp = reshape(inferedNetwork',NodesPerRowCol,NodesPerRowCol);
                                imagesc(inferedNetworkForDisp)
                                title(['inferedNetwork' ...
                                    ', Iter #' num2str(curIteration) ...
                                    ', CorrectLvl: ' num2str(correctLevel)]);

                                subplot(4,1,4)
                                imagesc(storedPatterns{inputPattern}-inferedNetworkForDisp)
                                title(['DiffFromAll[' num2str(inputPattern) ']' ...
                                    ', Iter #' num2str(curIteration) ...
                                    ', CorrectLvl: ' num2str(correctLevel)]);
                                
                                if (SAVE_PLOTS)
                                    saveas(h, ['Async' ...
                                        '_Pat' num2str(inputPattern) ...
                                        iif(DIVERT, ['_Divert' num2str(curDivert)], ...
                                            ['_Noise' num2str(curNoise)]) ...
                                        '~' datestr(now,'dd_mm_yy_HH_MM_SS_FFF') ...
                                        '.jpg']);
                                end
                            end
                        end %if PLOT_DIFF
                    
                        curIteration=curIteration+1;
                    end %while iterations
                end %if SYNC_UPDATE
            end %for t:curNumOfPatterns
            
            fprintf('Correct: %d, Iters: %d\n', ...
                numOfIterationsUntilInf(curNoiseLevel,curNumOfPatterns), ...
                numOfIterationsUntilCorrectInf(curNoiseLevel,curNumOfPatterns));
        end %for curNumOfPatterns
        curNoise = curNoise+NOISE_STEP;
        noisyNodesNumber= NodesPerPattern*curNoise;
    end
    
    %reset noise
    curNoise = INITIAL_NOISE_PERCENT;
    
    curDivert=curDivert+1;
end


%% NvsP Iterations
xspan=1:length(numOfIterationsUntilCorrectInf);
%yspan = ceil(numOfIterationsUntilCorrectInf./numOfCorrectInferences);
yspan = ceil(numOfIterationsUntilInf./numOfCorrectInferences);

plot(xspan,yspan, 'b');
xlabel('Patterns Stored');
ylabel('Patterns Infered');
title(['NotNaturalImages' ' - Average Iterations' ' - ' ...
    iif(SYNC_UPDATE, '', 'A') ...
    'Sync Update']);

t=toc;
disp(['Elapsed time is ' datestr(datenum(0,0,0,0,0,t),'HH:MM:SS') '.']);

beep;
