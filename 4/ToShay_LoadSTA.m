clc;
clear all;
close all;
figure;

for RepFlag=1:2 %1 Rep, 2 Non-Rep
    cell=1;

    load(['MatFiles\STA_Cell_' num2str(cell) ...
        '_RepFlag_' num2str(RepFlag) '.mat']);
   
    
    %normalize
    STA=accSTA/countSpikes;
    
    if RepFlag==1
        meanSTA_Rep = mean(STA);
        Normalized_STA_Rep=(STA-meanSTA_Rep);
        maxSTA_Rep=abs(max(Normalized_STA_Rep));
        Normalized_STA_Rep=Normalized_STA_Rep./maxSTA_Rep;
        accSTADS_Rep=(accSTADS_Rep-meanSTA_Rep)./maxSTA_Rep;
    else
        meanSTA_NonRep = mean(STA);
        Normalized_STA_NonRep=(STA-meanSTA_NonRep);
        maxSTA_NonRep=abs(max(Normalized_STA_NonRep));
        Normalized_STA_NonRep=Normalized_STA_NonRep./maxSTA_NonRep;
        accSTADS_NonRep=(accSTADS_NonRep-meanSTA_NonRep)./maxSTA_NonRep;
    end
    
    if RepFlag==1
        p1=plot([-SIM_TICKS_PER_SECOND:-1]/AP_TICKS_PER_SECOND,(Normalized_STA_Rep),'k');hold on      
    else
        p2=plot([-SIM_TICKS_PER_SECOND:-1]/AP_TICKS_PER_SECOND,(Normalized_STA_NonRep),'r');hold on
        xlabel('Time before AP (sec)');ylabel('STA');
        title(['Cell num. ' num2str(iNeuron)]);
        xlim([-SIM_TIME_FACTOR 0]);
        legend('StimulusRep','StimulusNonRep','location','northwest')      
    end  
    
end

%% Linear filter phase
TRAIL_IN_SECONDS = 100;
TRAIL_IN_AP_TICKS = TRAIL_IN_SECONDS*AP_TICKS_PER_SECOND;
TRAIL_IN_SIM_TICKS=TRAIL_IN_SECONDS*SIM_TICKS_PER_SECOND; %100 seconds (when we have 5000 per second)

%Normalizing by non-rep mean and max
NormalizedStimuliNonRep=(StimulusNonRep-meanSTA_NonRep)/(maxSTA_NonRep);

Filter = Normalized_STA_NonRep;
block=zeros(2,4000); %range is big enough size for remap, can be 4k

sumAbsSTA=sum(abs(Filter)); %convolution will be normalized by this, abs(STA window)

%r=0; %unused
EndLoop=floor(LightChange(end)/TRAIL_IN_SIM_TICKS)*TRAIL_IN_SIM_TICKS; %takes only full light times (removes inbetween ticks)
Filter=fliplr(Filter'); %flips the window (because it's not a conv, it's a sliding window/cross corr)

%iterate in steps of 100-1 seconds on stims
i = 0; %used for tracing
%jumping in trails (actually starting a second before each trail starts)
for t=LightChange(1):TRAIL_IN_SIM_TICKS-(1*SIM_TICKS_PER_SECOND):EndLoop 
    i=i+1;
    if (mod(i, 20)==0)
        fprintf('conving: %d/%d ...\n', round(t), EndLoop); %~
    end
        
    idxLightChange=find(LightChange<=t,1,'last'); %because we floored it, we look for the closest light change
    %we start right after
        
    %calc the light changes in trail, but in stim time scale
    %LIGHT_CHANGES_IN_TRAIL = TRAIL_IN_SIM_TICKS/10000*30;
    LIGHT_CHANGES_IN_TRAIL = (TRAIL_IN_SECONDS*SIM_TIME_FACTOR) ...
        * (1/STIMS_PER_SECOND);
    
    %get 100 seconds stims (appear every 1/30 sec) + 10 safety light changes (will be removed later)
    NormalizedStimuliInTrail=NormalizedStimuliNonRep(idxLightChange+1:idxLightChange+LIGHT_CHANGES_IN_TRAIL+10);
    %smear the light changes
    NormalizedStimuliInTrail=repmat(NormalizedStimuliInTrail',ceil(STIMS_PER_SECOND_IN_TICKS),1);
    NormalizedStimuliInTrail=NormalizedStimuliInTrail(:);
    
    %Removes every stim suffix added due to ceil, expect the 3rd stim
    %suffix, because we assume it falls on an interger, see attached photo1
    NormalizedStimuliInTrail(ceil(STIMS_PER_SECOND_IN_TICKS)*2:ceil(STIMS_PER_SECOND_IN_TICKS)*3:end)=[]; %removes every 3rd stim suffix, from the 2nd stim suffix
    NormalizedStimuliInTrail(ceil(STIMS_PER_SECOND_IN_TICKS):(ceil(STIMS_PER_SECOND_IN_TICKS)-1)+ceil(STIMS_PER_SECOND_IN_TICKS)*2:end)=[]; %removes every 3rd stim suffix, from the 1st stim suffix
    
    toRemove=t-LightChange(idxLightChange); %get only trail data
    NormalizedStimuliInTrail(1:round(toRemove))=[];
    NormalizedStimuliInTrail(TRAIL_IN_SIM_TICKS+1:end)=[];
    
    Conv=conv(Filter,NormalizedStimuliInTrail)/sumAbsSTA; %do the conv, and normalize by sum(abs(STA window))
    Conv(1:length(Filter))=[]; %throw a whole kernel window from start
    Conv(TRAIL_IN_SIM_TICKS-length(Filter)+1:end)=[]; %throw a whole kernel window from end
    StimXSta=Conv;
    
    %linear transform: remaping values to positive indexes
    X_ind=round((StimXSta+6)*10^2);
    
    %there are 2 ways to create rate:
    %1. using bins
    %2. using a gausian sliding window (current implementation)
    SpikeLen=1000;
    shift=length(Filter)/2; %because it's a low pass filter, we must shift by half and reduce half at the end
    %TODO: investigate, we might make simpler
    SpikeTime=find(TTNonRep(1,iNeuron).sp>=shift+t+(shift-1)-SpikeLen/2 & ...
        TTNonRep(1,iNeuron).sp<shift+t+TRAIL_IN_SIM_TICKS-1-shift+SpikeLen/2);
    SpikeVec=zeros(1,TRAIL_IN_SIM_TICKS-SIM_TICKS_PER_SECOND+SpikeLen);
    SpikeVec(round(TTNonRep(1,iNeuron).sp(SpikeTime) - ...
        (shift+t+(shift-1)-SpikeLen/2)+1 )) = 1;
    gaussWinFilter=gausswin(SpikeLen);
    SpikeRate=conv(gaussWinFilter,SpikeVec)/sum(gaussWinFilter);
    SpikeRate(1:SpikeLen)=[];
    SpikeRate(TRAIL_IN_SIM_TICKS-length(Filter)+SpikeLen-SpikeLen+1:end)=[];
    
    %append rate for the corresponding X vals
    block(1,X_ind')= block(1,X_ind')+SpikeRate;
    block(2,X_ind')= block(2,X_ind')+1;
end

%removeFirst=find(sum(block)>0,1,'first'); unused
%removeLast=find(sum(block)>0,1,'last'); %unused

Y=block(1,:)./(block(2,:)+.000001); %calc the rate using division of counters
%add .000001 so we won't have NaNs (division by zero = NaN)

X=[1:size(Y,2)]/10^2-6; %linear transform: remaping values to positive indexes

remove=find(block(2,:)<5);%Here take off noise 
Y(remove)=[];
X(remove)=[];

if iNeuron==2 || iNeuron==3 %TODO: check if in use
    remove=find(X<-1); 
    Y(remove)=[];
    X(remove)=[];
end

[func, gof] =FitLNmodel(X,Y,iNeuron);%Y units:1/sample=10,000Hz
ylabel('FR [Hz*10,000]');xlabel('Stim*STA');hold on
legend( 'From RD', 'Fit' );%raw data

%% using the stimulus rep
NormalizedStimuliRep=(StimulusRep-meanSTA_NonRep)/(maxSTA_NonRep);%normalized by non rep params

%OLD:TRAIL_IN_SIM_TICKS=10*10000;%500000;
TRAIL_IN_SIM_TICKS = TRAIL_IN_SECONDS*SIM_TICKS_PER_SECOND; %100 seconds (when we have 5000 per second)

EndLoop=floor(LightChangeRep(end)/TRAIL_IN_SIM_TICKS)*TRAIL_IN_SIM_TICKS; %takes only full light times (removes inbetween ticks)
t=TRAIL_IN_SECONDS*AP_TICKS_PER_SECOND; %start after the first trail
idxLightChange=find(LightChangeRep<=t,1,'last');

LIGHT_CHANGES_IN_TRAIL = (TRAIL_IN_SECONDS*SIM_TIME_FACTOR) ...
        * (1/STIMS_PER_SECOND);

%get 100 seconds stims (appear every 1/30 sec) + 10 safety light changes (will be removed later)
NormalizedStimuliInTrail=NormalizedStimuliRep(idxLightChange+1:idxLightChange+LIGHT_CHANGES_IN_TRAIL+10);

%{
OLD
%get 100 seconds stims (appear every 1/30 sec) + 10 safety light changes (will be removed later)
NormalizedStimuliInTrail=NormalizedStimuliRep(idxLightChange+1:idxLightChange+round(TRAIL_IN_SIM_TICKS/10000*30)+10);%all colors in the current 5000 window
%}

%smear the light changes
NormalizedStimuliInTrail=repmat(NormalizedStimuliInTrail',ceil(STIMS_PER_SECOND_IN_TICKS),1);
NormalizedStimuliInTrail=NormalizedStimuliInTrail(:);
    
%Removes every stim suffix added due to ceil, expect the 3rd stim
%suffix, because we assume it falls on an interger, see attached photo1
NormalizedStimuliInTrail(ceil(STIMS_PER_SECOND_IN_TICKS)*2:ceil(STIMS_PER_SECOND_IN_TICKS)*3:end)=[]; %removes every 3rd stim suffix, from the 2nd stim suffix
NormalizedStimuliInTrail(ceil(STIMS_PER_SECOND_IN_TICKS):(ceil(STIMS_PER_SECOND_IN_TICKS)-1)+ceil(STIMS_PER_SECOND_IN_TICKS)*2:end)=[]; %removes every 3rd stim suffix, from the 1st stim suffix
    
toRemove=t-LightChangeRep(idxLightChange); %get only trail data
NormalizedStimuliInTrail(1:round(toRemove))=[];
NormalizedStimuliInTrail(TRAIL_IN_SIM_TICKS+1:end)=[];

Conv=conv(Filter,NormalizedStimuliInTrail)/sumAbsSTA; %do the conv, and normalize by sum(abs(STA window))
Conv(1:length(Filter))=[]; %throw a whole kernel window from start
Conv(TRAIL_IN_SIM_TICKS-length(Filter)+1:end)=[]; %throw a whole kernel window from end
StimXSta=Conv;
%Execute the Non linear fit
SimulatedFR=func(StimXSta);


StimTimeRepFixed=StimTimeRep(1:2:end); %fix stimTimeRep to be 200 sec diff
TT=TTRep;
APTimesCrossTrails=[];
countTrails=1; %counts the trails

REP_TRAIL_IN_SECONDS = 200;
REP_TRAIL_IN_AP_TICKS = REP_TRAIL_IN_SECONDS*AP_TICKS_PER_SECOND;

for trialOnset=StimTimeRepFixed(1:end-1)';
    %locate trail
    firstAP=find(TT(1,iNeuron).sp>trialOnset,1,'first');
    lastAP=find(TT(1,iNeuron).sp<=trialOnset+REP_TRAIL_IN_AP_TICKS,1,'last'); 
    %set all its aps times relative to the stimTime
    APTimesCrossTrails(end+1:end+lastAP-firstAP+1)=TT(1,iNeuron).sp(firstAP:lastAP)-trialOnset;
    countTrails=countTrails+1; 
end

%% plot

TICKS_PER_BIN = 100;
nbins=REP_TRAIL_IN_AP_TICKS/TICKS_PER_BIN; %100 ticks each bin
SpikeLen=10; %how many spikes we consider 1 spike (smoothing), ORIGINAL=10
figure;

[binCounts,xCenters]=hist(APTimesCrossTrails/AP_TICKS_PER_SECOND,nbins);
hold on;

gaussWinFilter=gausswin(SpikeLen);
actualRate=conv(gaussWinFilter,binCounts/countTrails)/sum(gaussWinFilter);%units:1/sample
%lowpassfilter: start after half a filter window and end before half
actualRate = actualRate(SpikeLen/2:end-SpikeLen/2);


%plot(xCenters,binCounts/countTrails); %actual spike train
plot(xCenters,actualRate,'r'); 

ylabel('FR [10,000*Hz]');
xlim([t/AP_TICKS_PER_SECOND t/AP_TICKS_PER_SECOND+10]);

hold on
shift=length(Filter)/2; %because it's a low pass filter, we must shift by half and reduce half at the end
estRate = SimulatedFR*100;
plot((0+(t+shift:t+TRAIL_IN_SIM_TICKS-(shift+1)))/AP_TICKS_PER_SECOND,estRate,'k');
xlabel('Time [sec]');
title(['Cell num. ' num2str(iNeuron)]);
legend('PSTH','Simulated FR')


%% STC
figure;
A=downsample(Normalized_STA_Rep,DOWN_SAMPLE_EVERY_NTH);%200 instead of 5000
N_spikes=size(accSTADS_NonRep,1);
STC=1/(N_spikes-1)*((accSTADS_NonRep'-repmat(A,1,N_spikes))*...
    (accSTADS_NonRep'-repmat(A,1,N_spikes))');
time=[1:200]/400;
colormap(gray);
surf(repmat(time',1,200),repmat(time,200,1),STC);

grid off;shading flat
set(gca,'view',[0 90])
colorbar
title(['Cell num. ' num2str(iNeuron)]);


figure;

plot(eig(STC),'.')
xlim([0 200])
xlabel('Eigenvalue number');ylabel('Eigenvalue (Variance)');
title(['Cell num. ' num2str(iNeuron)]);

[V,D]=eig(STC);

figure
plot(V(:,200));
title(['Cell num. ' num2str(iNeuron)]);