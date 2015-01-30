load Data.mat
%convert to 1/30000 sec
time = 3;

numberOfNeurons=4;

noOfstimuliPerRun = 6000;
timeBins =30000;
numberOfruns=200;
staLength=12000; %12000 equals 400 m.s.


staStack = zeros(staLength,numberOfNeurons);

figure;


stimtimeNonRepReduced = StimTimeNonRep(1:2:end)*time;


noOfReapeats=10;
StimulusRepAfterConversion = StimulusNonRep(:,ones(noOfReapeats,1))';
StimulusRepAfterConversion=StimulusRepAfterConversion(:);


s = noOfstimuliPerRun;
tempStim = zeros(s,1);
for ii=1:numberOfNeurons
    subplot(numberOfNeurons+1,1,ii);
    
    hold on
    staSum=zeros(staLength,1);
    staCounter=zeros(staLength,1);
    
    n = TTNonRep(1,ii);
    ap=n.sp*time;
    ap(ap < stimtimeNonRepReduced(1,1)+staLength) = []; %throw anything which doesn't have a full stim
    counter=0;
    for i=1:length(ap)
        
        
        currentAp = ap(i);
        
        staBegin = currentAp-staLength;
        
        ind = find(stimtimeNonRepReduced(:,1)<currentAp,1,'last');
        if(staBegin<stimtimeNonRepReduced(ind,1))
         
            
        else
            
            tempStim= [stimtimeNonRepReduced(ind,1):1000:(stimtimeNonRepReduced(ind,1)+(timeBins*numberOfruns))-1]';
            
            currentWindow = currentAp-staLength:currentAp-1;
            
            OUTIND = find(ismember(tempStim,currentWindow));
            OUTIND2 = find(ismember(currentWindow,tempStim));
            if (ind>1)
                OUTIND=OUTIND+(6000*(ind-1));
            end
            %too many stimuli iterations
            dd =find(OUTIND>324000);
            if(isempty(dd))
                if(length(OUTIND)==11)
                    temp=zeros(12,1);
                    temp(2:12,1)=OUTIND;
                    temp(1,1)=temp(2,1)-1;
                    OUTIND=temp;
                    
                end
                s=StimulusNonRep(OUTIND);
                s = s(:,ones(1000,1))';
                s=s(:);
                if(length(s)~=12000)
                    lkj=8;
                end
                staSum = staSum+s;
                
            end
        end
        
        
    end
    
    staSum = staSum ./ length(ap);
    staStack(:,ii) = staSum;
%     plot([-staLength+1:0],staSum,'LineWidth',2);
%     xlabel('Time t in 1/30000 sec');
%     xlim([-staLength 0]);
%     ylim([0.35 0.55]);
%     ylabel('STA');
  stemp=staSum(1:1000:12000);
    %figure;
    plot((1:12)/30, stemp , '--');
    title(['sta neuron ' num2str(ii)]);
    xlabel('Time (1/1000 sec)')
    
%     
%     
%     %calc nonlinear filter (by ratio)
%     sTemp = staSum(1:1000:end);
%     f1=filter(sTemp,1,StimulusNonRep);
%     hist(f1(13:end))
%     
%     apPreSta = zeros(length(ap),1);
%     for i=1:length(ap)
%         
%         
%         currentAp = ap(i);
%         
%         staBegin = currentAp-staLength;
%         
%         ind = find(stimtimeNonRepReduced(:,1)<currentAp,1,'last');
%         if(staBegin<stimtimeNonRepReduced(ind,1))
%             
%         else
%             
%             tempStim= [stimtimeNonRepReduced(ind,1):1000:(stimtimeNonRepReduced(ind,1)+(timeBins*numberOfruns))-1]';
%             tempStimAfterTimeConversion = tempStim(:,ones(5,1))';
%             tempStimAfterTimeConversion=tempStimAfterTimeConversion(:);
%             currentWindow = currentAp-staLength:currentAp-1;
%             
%             OUTIND = find(ismember(tempStim,currentWindow));
%             % OUTIND2 = find(ismember(currentWindow,tempStim));
%             if (ind>1)
%                 OUTIND=OUTIND+(6000*(ind-1));
%             end
%             %too many stimuli iterations
%             dd =find(OUTIND>324000);
%             if(isempty(dd))
%                 if(length(OUTIND)==11)
%                     temp=zeros(12,1);
%                     temp(2:12,1)=OUTIND;
%                     temp(1,1)=temp(2,1)-1;
%                     OUTIND=temp;
%                     
%                 end
%                 s=StimulusNonRep(OUTIND);
%                 
%                 
%                 apPreSta(i) = dot(s,sTemp);
%                 
%                 
%             end
%             
%         end
%         
%     end
%     
%     
    
end
raw = f1(13:end);
[n1,h1]  = hist(raw,12);

[n2,h2] = hist(apPreSta,12);
ratio = h2./h1;
plot(ratio);
subplot(numberOfNeurons+1,1,numberOfNeurons+1);
bar(0:diff,1,'r')
% plot([0:diff diff:(diff+diff)],[zeros(1000,1) ones(1000,1)],'r')
axis([-1 100 0 1])
xlabel('Light switch')
ylabel('diffrent light');
s=1;