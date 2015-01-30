load Data.mat
measure=30000;
%convert to 1/30000 sec
time = 3;
noOfReapeats = 1000;
numberOfNeurons=4;
diff=600000; %1200 equals 40 m.s.
stimWindow = 600000;
noOfstimuliPerRun = 6000;
timeBins =60000;

bin=12000;

noBin =stimWindow/bin;

hrtz = 1;


StimTimeRepEvery200 = StimTimeRep(1:2:end)*time;



%{
StimulusRepAfterTimeConversion = StimulusRep(:,ones(noOfReapeats,1))';
StimulusRepAfterTimeConversion=StimulusRepAfterTimeConversion(:);
%}

psthSum = zeros(noBin,numberOfNeurons);
for ii=1:numberOfNeurons
    subplot(numberOfNeurons+1,1,ii);
    
    hold on
    
    psth = zeros(noBin,1); %init psth
    for i=1:length(stimtimeRepReduced)
        subStimuli = 0;
        
        n = TTRep(1,ii);
        ap=n.sp*time;
        
        
        currentOff = stimtimeRepReduced(i);
        %current chunk of 200 secs
        
        %on = currentOff-diff;
        nextOn= currentOff+diff;
        %look until next chunk/repetition
        
        
        IDX = uint32(1:size(ap,1));
        ind = IDX(ap(:,1) >= currentOff & ap(:,1) < nextOn);
        currentAP=ap(ind);
        %get all APs in the chunk
        
      
       scatter((currentAP-currentOff)/3, ones(length(currentAP),1)*i,1,'b');
      %get back to 10000 scale
        
        
        
        
        
        
        
        
    end
    ylim([0 22])
    ylabel('Trial number')
    xlabel('Time (1/10000 sec)')
    
end
psthSum = psthSum
subplot(numberOfNeurons+1,1,numberOfNeurons+1);
bar(0:diff:diff,1,'r')
% plot([0:diff diff:(diff+diff)],[zeros(1000,1) ones(1000,1)],'r')
%axis([-1 stimWindow 0 1])
xlabel('Light switch')
ylabel('diffrent light');
hjs=1;