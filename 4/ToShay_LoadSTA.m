clc;
clear all;
close all;
figure;%subplot(4,1,cell);

for RepFlag=1:2 %1 Rep, 2 Non-Rep
    cell=1;

    load(['MatFiles\STA_Cell_' num2str(cell) ...
        '_RepFlag_' num2str(RepFlag) '.mat']);
   
    
    if RepFlag==1
        staMeanRep=[]; %mean over APs
        meansta=[]; %mean over all
        maxsta=[]; %max abs
        NstaMeanRep=[]; %normalized by max, will be ploted and used as filter
        
        staMeanRep(:,cell)=staSum/(i-firstAP+1);
        meansta=mean(staMeanRep(:,cell));
        staZ(:,cell)=staMeanRep(:,cell)-meansta;
        maxstaZ=abs(max(staZ(:,cell)));
        NstaMeanRep(:,cell)=(staZ(:,cell))/(maxstaZ);
        sMatRep=(sMatRep-meansta)./(maxstaZ); %normalize
        
    else
        staMeanNonRep=[]; %mean over APs
        meanstaNon=[]; %mean over all
        maxstaNon=[]; %max abs
        NstaMeanNonRep=[]; %normalized by max, will be ploted and used as filter
     
        staMeanNonRep(:,cell)=staSum/(i-firstAP+1);
        meanstaNon=mean(staMeanNonRep(:,cell));
        staZNon(:,cell)=staMeanNonRep(:,cell)-meanstaNon;
        maxstaZNom=abs(max(staZNon(:,cell)));
        NstaMeanNonRep(:,cell)=(staZNon(:,cell))/(maxstaZNom);
        sMatNonRep=(sMatNonRep-meanstaNon)/(maxstaZNom); %normalize
    end
    
end

for RepFlag=1:2 %1 Rep, 2 Non-Rep
    if RepFlag==1
            p1=plot([-time*10000:-1-a]/10000,(NstaMeanRep(:,cell)),'k');hold on      
        else
            p2=plot([-time*10000:-1-a]/10000,(NstaMeanNonRep(:,cell)),'r');hold on
            xlabel('Time before AP (sec)');ylabel('STA');
            title(['Cell num. ' num2str(cell)]);
            xlim([-time 0]);
            legend([p1 p2],'StimulusRep','StimulusNonRep','location','northwest')      
    end  
end

%%        
N_stim=(StimulusNonRep-meanstaNon)/(maxstaZNom);%normalized stims
N_stimRep=(StimulusRep-meanstaNon)/(maxstaZNom);%normalized stims
window= NstaMeanNonRep(:,cell);
block=zeros(2,40000);
%blockR=zeros(50,40000); %unused
total1=ones(1,5000)*(abs(window)); %convolution will be normalized by this, abs(STA window)
%X_ind=[];FR=[];xSTA=[];nStim=[]; %useless
ConvLen=500000; %100 seconds (when we have 500 per second)
%r=0; %unused
EndLoop=floor(LightChange(end)/ConvLen)*ConvLen; %takes only full light times (removes inbetween ticks)
window=fliplr(window'); %flips the window (because it's not a conv, it's a sliding window/cross corr)

%iterate in steps of 100-1 seconds on stims
i = 0; %used for tracing
convsRange = LightChange(1):ConvLen-5000:EndLoop;
for t=convsRange
    i=i+1;
    if (mod(i, 10)==0)
        fprintf('conving: %d/%d ...\n', i, length(convsRange));
    end
    
    %r=r+1; %unused
    ind=find(LightChange<=t,1,'last'); %because we floored it, we look for the closest light change
    %we start right after
    
    %toSample=N_stim(ind+1:ind+round(ConvLen/10000*30)+10);%round in useless
    toSample=N_stim(ind+1:ind+ConvLen/10000*30+10);%all collors in the current 5000 window
    %get 100 seconds stims (appear every 1/30 sec) + 10(/5000) safety ticks
    
    N_stim_sampled=repmat(toSample',ceil(10000/30),1);
    N_stim_sampled=N_stim_sampled(:);
    
    %Removes every stim suffix added due to ceil, expect the 3rd stim
    %suffix, because we assume it falls on an interger, see attached photo1
    N_stim_sampled(668:334*3:end)=[]; %removes every 3rd stim suffix, from the 2nd stim suffix
    N_stim_sampled(334:333+334*2:end)=[]; %removes every 3rd stim suffix, from the 1st stim suffix
    
    toRemove=t-LightChange(ind); %get only trail data
    N_stim_sampled(1:round(toRemove))=[];  %get only light change after t
    N_stim_sampled(ConvLen+1:end)=[];
    
    Conv=conv(window,N_stim_sampled)/total1; %do the conv, and normalize by abs(STA window)
    Conv(1:5000)=[]; %throw a whole STA window from start
    Conv(ConvLen-5000+1:end)=[]; %throw a whole STA window from end
    StimXSta=Conv;
    
    X_ind=round((StimXSta+6)*10^2);%bin scale
    
    %there are 2 ways to create rate:
    %1. using bins
    %2. using a gausian sliding window (current implementation)
    SpikeLen=1000;shift=2500;
    SpikeTime=find(TTNonRep(1,cell).sp>=shift+t+2499-SpikeLen/2 & TTNonRep(1,cell).sp<shift+t+ConvLen-1-2500+SpikeLen/2);
    SpikeVec=zeros(1,ConvLen-5000+SpikeLen);
    SpikeVec(round(TTNonRep(1,cell).sp(SpikeTime) -(shift+t+2499-SpikeLen/2)+1 ))=1;
    a=gausswin(SpikeLen);
    SpikeRate=conv(a,SpikeVec)/sum(a);
    SpikeRate(1:SpikeLen)=[];
    SpikeRate(ConvLen-5000+SpikeLen-SpikeLen+1:end)=[];
    
    %append rate for the corresponding X vals
    block(1,X_ind')= block(1,X_ind')+SpikeRate;
    block(2,X_ind')= block(2,X_ind')+1;
end

removeFirst=find(sum(block)>0,1,'first');
removeLast=find(sum(block)>0,1,'last');

Y=block(1,:)./(block(2,:)+.000001); %calc the rate using division of counters
X=[1:size(Y,2)]/10^2-6;

remove=find(block(2,:)<5);%Here take off noise 
Y(remove)=[];
X(remove)=[];

if cell==2 || cell==3 %TODO: check if in use
    remove=find(X<-1); 
    Y(remove)=[];
    X(remove)=[];
end

[func, gof] =FitLNmodel(X,Y,cell);%Y units:1/sample=10,000Hz
ylabel('FR [Hz*10,000]');xlabel('Stim*STA');hold on
legend( 'From RD', 'Fit' );%raw data

%% using the stimulus rep

N_stim=(StimulusRep-meanstaNon)/(maxstaZNom);%normalized colors
window= NstaMeanNonRep(:,cell);
total1=ones(1,5000)*(abs(window));

%X_ind=[];FR=[];xSTA=[];nStim=[]; %unused
ConvLen=10*10000;%500000;
EndLoop=floor(LightChangeRep(end)/ConvLen)*ConvLen;
vec=[];
window=fliplr(window');
t=100*10000;

ind=find(LightChangeRep<=t,1,'last');

toSample=N_stim(ind+1:ind+round(ConvLen/10000*30)+10);%all colors in the current 5000 window
N_stim_sampled=repmat(toSample',ceil(10000/30),1);
N_stim_sampled=N_stim_sampled(:);
N_stim_sampled(668:334*3:end)=[];
N_stim_sampled(334:333+334*2:end)=[];
toRemove=t-LightChangeRep(ind);
N_stim_sampled(1:round(toRemove))=[];
N_stim_sampled(ConvLen+1:end)=[];

Conv=conv(window,N_stim_sampled)/total1;
Conv(1:5000)=[];
Conv(ConvLen-5000+1:end)=[];
StimXSta=Conv;
SimulatedFR=func(StimXSta);

StimTime1=StimTimeRep(1:2:end);
TT=TTRep;
AP=[];
r=1; %counts the stims
row=[]; %TODO: unused
for lo=StimTime1(1:end-1)';
    %locate trail
    First=find(TT(1,cell).sp>lo,1,'first');
    Last=find(TT(1,cell).sp<=lo+2000000,1,'last'); 
    %set all its aps relative to the stimTime
    AP(end+1:end+Last-First+1)=TT(1,cell).sp(First:Last)-lo;
    row(end+1:end+Last-First+1)=r; %TODO: unused
    r=r+1; 
end

nbin=20000;
SpikeLen=10;
figure;

[y x]=hist(AP/10000,nbin);hold on

yGauss=conv(gausswin(SpikeLen),y/r)/sum(gausswin(SpikeLen));%units:1/sample
plot(x,yGauss(SpikeLen/2:end-SpikeLen/2),'r')
ylabel('FR [10,000*Hz]');
xlim([t/10000 t/10000+10])

hold on
plot((0+[t+2500:t+ConvLen-2501])/10000,SimulatedFR*100,'k');
xlabel('Time [sec]');
title(['Cell num. ' num2str(cell)]);
legend('PSTH','Simulated FR')

%% STC

A=downsample(NstaMeanRep(:,cell),25);%200 instead of 5000
N_spikes=size(sMatNonRep,1);
STC=1/(N_spikes-1)*((sMatNonRep'-repmat(A,1,N_spikes))*...
    (sMatNonRep'-repmat(A,1,N_spikes))');
time=[1:200]/400;
colormap(gray);
surf(repmat(time',1,200),repmat(time,200,1),STC);

grid off;shading flat
set(gca,'view',[0 90])
colorbar
title(['Cell num. ' num2str(cell)]);
figure;

plot(eig(STC),'.')
xlim([0 200])
xlabel('Eigenvalue number');ylabel('Eigenvalue (Variance)');
title(['Cell num. ' num2str(cell)]);

[V,D]=eig(STC);

figure
plot(V(:,200));
title(['Cell num. ' num2str(cell)]);

load gong 
sound(y,Fs)