load Data.mat
measure=30000;
%convert to 1/30000 sec
time = 3;
noOfReapeats = 1000;
numberOfNeurons=4;
diff=6000000; %1200 equals 40 m.s.
stimWindow = 12000000;
noOfstimuliPerRun = 6000;
timeBins =60000;

bin=100000;

noBin =stimWindow/bin;

hrtz = 1;


stimtimeRepReduced = StimTimeRep(1:2:end)*time;




StimulusRepAfterTimeConversion = StimulusRep(:,ones(noOfReapeats,1))';
StimulusRepAfterTimeConversion=StimulusRepAfterTimeConversion(:);


psthSum = zeros(noBin/2,numberOfNeurons);
for ii=1:numberOfNeurons
    subplot(numberOfNeurons+1,1,ii);
    
    hold on
    
    psth = zeros(noBin,1); %init psth
    for i=1:length(stimtimeRepReduced)
        subStimuli = 0;
        
        n = TTRep(1,ii);
        ap=n.sp*time;
        
        
        currentOff = stimtimeRepReduced(i);
        
        on = currentOff-diff;
        nextOn= currentOff+diff;
        
        
        IDX = uint32(1:size(ap,1));
        ind = IDX(ap(:,1) >= on & ap(:,1) < nextOn);
        
        
        fromBin=0;
        apBin= ap(ind)-on;
        for b=1:noBin
            toBin = fromBin+bin;
            IDX = uint32(1:size(apBin,1));
            indBin = IDX(apBin(:,1) >= fromBin & apBin(:,1) < toBin);
            psth(b,1) = psth(b,1)+length(apBin(indBin));
            fromBin = toBin;
        end
        
        
        
        
        
        
        
        
    end
    
    if(hrtz==1)
        psthSum(:,ii) = psth(61:120);
        
        bar(0:bin:(stimWindow-bin),(psth/(noOfstimuliPerRun*length(stimtimeRepReduced)*bin))*diff)
        
        ylabel('Rate(Hz)');
        axis([-1 stimWindow 0 max((psth/(noOfstimuliPerRun*length(stimtimeRepReduced)*bin))*diff)])
    else
        bar(0:bin:(stimWindow-bin),psth)
        
        axis([-1 stimWindow 0 max(psth)])
    end
    xlabel(['Time (12000000= 400 sec) , Bin size in sec: ' num2str(bin/measure) ] )
    
    
end
psthSum = psthSum
subplot(numberOfNeurons+1,1,numberOfNeurons+1);
bar(0:diff:diff,1,'r')
% plot([0:diff diff:(diff+diff)],[zeros(1000,1) ones(1000,1)],'r')
%axis([-1 stimWindow 0 1])
xlabel('Light switch')
ylabel('diffrent light');
hjs=1;