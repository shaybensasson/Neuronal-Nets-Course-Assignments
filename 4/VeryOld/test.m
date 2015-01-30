close all;
DATA_POINTS=50;
figure;
timeOfAPs=TTNonRep(1).sp*TIME_FACTOR;
X=StimTimeNonRepFixed(53):2000:StimTimeNonRepFixed(54);
Y=StimulusNonRep(53*2000+1:54*2000+1)';
X = X(1:DATA_POINTS);
Y = Y(1:DATA_POINTS);
scatter(X, Y);
hold on;

%set(gca,'XTickLabel',X);
%pause;
medY = (min(Y)+max(Y))/2;
minY = min(Y);

start = min(X);
sof = max(X);

indexes = find(timeOfAPs>=start & timeOfAPs<=sof,DATA_POINTS,'first');
length(indexes)
X2 = timeOfAPs(indexes);
Y2 = ones(length(indexes),1)*medY;
scatter(X2, Y2, 40, 'r');
iAPTime = indexes(end);

curAPTime = timeOfAPs(iAPTime);
        
        %we look back STA_WINDOW time
        staBegin = curAPTime;
        
        %indexes = find(StimTimeNonRepFixed(:)<curAPTime,1,'last');
        
        
        %lets find the last stimulus that occured before the AP
        iClosestStimTime = find(StimTimeNonRepFixed(:)<curAPTime,1,'last');
        %{
        X3 = StimTimeNonRepFixed(iClosestStimTime):2000:curAPTime;
        Y3 = ones(length(X3),1)*minY;
        scatter(X3, Y3, 40,'y', 'fill');
        %}
        
        if(staBegin<StimTimeNonRepFixed(iClosestStimTime))
            %the closest stimulus we found is inside the window, 
            %we shall discard the ap (too fast for our STA window).
            %We'd like that STA will begin inside the window, 
            %so we can find intersection stims and the ap
            %{
            fprintf( ...
                        ['staBegin<=StimTimeNonRepFixed(iClosestStimTime)' ...
                        ' N(%d),APTime(%d)\n'], ...
                           iNeuron, iAPTime);
            %}
        else %STA begins inside stim time non-repetition window
            
            %define a window backwards in time STA_WINDOW long
            %we take -1, because of causality - stim causes ap
            currentWindow = curAPTime-STA_WINDOW:curAPTime-1;
            
            %get the stim non-repetition window starting from the closest stim time 
            % (jumping in STIMULUS_INTERVAL, occurs every 1/30 sec)
            stimuliTimes = ... 
                (StimTimeNonRepFixed(iClosestStimTime):STIMULUS_INTERVAL: ...
                (StimTimeNonRepFixed(iClosestStimTime)+NON_REPEAT_WINDOW)-1)';
            
            indexesOfStimTimesAtSTA = find(ismember(stimuliTimes,currentWindow));
            
            X3 = stimuliTimes(indexesOfStimTimesAtSTA);
            Y3 = ones(length(X3),1)*minY;
            scatter(X3, Y3, 40,'y', 'fill');
            
            %Shift the indexes by the iClosestStimTime index (fit index to
            % the approprite window)
            indexesOfStimTimesAtSTA = ...
                indexesOfStimTimesAtSTA + ...
                (LIGHT_LEVELS_BETWEEN_STIM_TIMES*(iClosestStimTime-1));
            
            %too many stimuli iterations
            %TODO: 324000/6000=54
            dd =find(indexesOfStimTimesAtSTA>324000);
            if(isempty(dd))
                countSimuli = length(indexesOfStimTimesAtSTA);
                if (countSimuli<EXPECTED_STIMULI_COUNT)
                        fprintf( ...
                            ['NOTE: stimuliAtSTA != %d: %d.' ...
                            ' N(%d),APTime(%d)\n'], ...
                               EXPECTED_STIMULI_COUNT, ...
                               countSimuli, iNeuron, iAPTime);
                end

                if(countSimuli==(EXPECTED_STIMULI_COUNT-1)) 
                    %we expect EXPECTED_STIMULI_COUNT, but we substract 1
                    % at the end of 'currentWindow', due to causality.
                    %ex. so instead of 12000/1000=12 will have 11 stimuli
                    % to solve this we add an earlier datapoint (index)
                    temp=zeros(EXPECTED_STIMULI_COUNT,1);
                    temp(2:EXPECTED_STIMULI_COUNT)=indexesOfStimTimesAtSTA;
                    temp(1)=temp(2)-1;
                    indexesOfStimTimesAtSTA=temp;
                    
                    %TODO: does it happen here also?
                end

                %get the stim values of the STA
                stimuliValuesAtSTA=StimulusNonRep(indexesOfStimTimesAtSTA);

                %{
                if(length(stimuliValuesAtSTA)~=EXPECTED_STIMULI_COUNT)
                    fprintf( ...
                        ['WARNING: length(stimuliValuesAtSTA)<%d>~=EXPECTED_STIMULI_COUNT<%d>.' ...
                        ' N(%d),APTime(%d)\n'], ...
                        length(stimuliValuesAtSTA), EXPECTED_STIMULI_COUNT, iNeuron, iAPTime);
                end
                %}

                %accumulate all values in STA window
                staSum = staSum+stimuliValuesAtSTA;
            end
        
        end
