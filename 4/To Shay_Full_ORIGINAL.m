%%
%STA 
staMeanRep=[];
staMeanNonRep=[];

cell=1
figure(5);%subplot(4,1,cell);
time=.5;%sec
a=0;
sMatRep=[];
sMatNonRep=[];

for RepFlag=1:2
    staSum=zeros(time*10000-a,1);
    StimTime=StimTimeRep;
    Stimulus=repmat(StimulusRep,length(StimTime(1:2:end)),1);%light intensity
    TT=TTRep;
    if RepFlag==2 %after calculating StimulusRep, now StimulusNonRep
        Stimulus=StimulusNonRep;%light intensity
        StimTime=StimTimeNonRep;
        TT=TTNonRep;
    end
    
    v1=length(StimTime);
    v2=[1:3000]*10000/30;%%%in the language of 10000 the colors
    v3=repmat(v2,v1,1); 
    LightChange=repmat(StimTime,1,3000);
    LightChange=v3'+LightChange';
    
    %changing the matrix to one long column
    LightChange=LightChange(:);
    LightChange(length(Stimulus)+1:end)=[];
    if RepFlag==1
        LightChangeRep=LightChange;
    end
    
    firstAP=find(TT(1,cell).sp>LightChange(1)+time*10000,1,'first');
    lastAP=find(TT(1,cell).sp<LightChange(end),1,'last');
    spike=0;
    for i=firstAP:lastAP%loop on ind of (TT) action potential
        ind=find(LightChange<TT(1,cell).sp(i)-time*10000,1,'last');
        ind2=find(LightChange<TT(1,cell).sp(i),1,'last');
        STA=Stimulus(ind+1:ind2+3);

        sta=repmat(STA',ceil(10000/30),1);
        sta=sta(:);
        sta(668:334*3:end)=[];
        sta(334:333+334*2:end)=[];
        toRemoveInit=round(TT(1,cell).sp(i)-time*10000-LightChange(ind));
        sta(1:toRemoveInit)=[];
        sta(time*10000+1-a:end)=[];
        spike=spike+1;
        DSsta = downsample(sta,25);
        
        if RepFlag==1
            sMatRep(spike,:)=DSsta;
        else
            sMatNonRep(spike,:)=DSsta;
        end
        staSum=staSum+sta;
    end

    if RepFlag==1
        staMeanRep=[];
        meansta=[];
        maxsta=[];
        NstaMeanRep=[];
        
        staMeanRep(:,cell)=staSum/(i-firstAP+1);
        meansta=mean(staMeanRep(:,cell));
        staZ(:,cell)=staMeanRep(:,cell)-meansta;
        maxstaZ=abs(max(staZ(:,cell)));
        NstaMeanRep(:,cell)=(staZ(:,cell))/(maxstaZ);
        sMatRep=(sMatRep-meansta)./(maxstaZ);
        
    else
        staMeanNonRep=[];
        meanstaNon=[];
        maxstaNon=[];
        NstaMeanNonRep=[];
     
        staMeanNonRep(:,cell)=staSum/(i-firstAP+1);
        meanstaNon=mean(staMeanNonRep(:,cell));
        staZNon(:,cell)=staMeanNonRep(:,cell)-meanstaNon;
        maxstaZNom=abs(max(staZNon(:,cell)));
        NstaMeanNonRep(:,cell)=(staZNon(:,cell))/(maxstaZNom);
        sMatNonRep=(sMatNonRep-meanstaNon)/(maxstaZNom);    
    end
    
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
N_stim=(StimulusNonRep-meanstaNon)/(maxstaZNom);%normalized colors
N_stimRep=(StimulusRep-meanstaNon)/(maxstaZNom);%normalized colors
window= NstaMeanNonRep(:,cell);
block=zeros(2,40000);
blockR=zeros(50,40000);
total1=ones(1,5000)*(abs(window));
X_ind=[];FR=[];xSTA=[];nStim=[];
flag=1;ConvLen=500000;r=0;
EndLoop=floor(LightChange(end)/ConvLen)*ConvLen;
window=fliplr(window');
for t=LightChange(1):ConvLen-5000:EndLoop%
    r=r+1;
    ind=find(LightChange<=t,1,'last');
    
    toSample=N_stim(ind+1:ind+round(ConvLen/10000*30)+10);%all collors in the current 5000 window
    N_stim_sampled=repmat(toSample',ceil(10000/30),1);
    N_stim_sampled=N_stim_sampled(:);
    N_stim_sampled(668:334*3:end)=[];
    N_stim_sampled(334:333+334*2:end)=[];
    toRemove=t-LightChange(ind);
    N_stim_sampled(1:round(toRemove))=[];
    N_stim_sampled(ConvLen+1:end)=[];
    
    Conv=conv(window,N_stim_sampled)/total1;
    Conv(1:5000)=[];
    Conv(ConvLen-5000+1:end)=[];
    StimXSta=Conv;
    
    X_ind=round((StimXSta+6)*10^2);%bin scale
    
    SpikeLen=1000;shift=2500;
    SpikeTime=find(TTNonRep(1,cell).sp>=shift+t+2499-SpikeLen/2 & TTNonRep(1,cell).sp<shift+t+ConvLen-1-2500+SpikeLen/2);
    SpikeVec=zeros(1,ConvLen-5000+SpikeLen);
    SpikeVec(round(TTNonRep(1,cell).sp(SpikeTime) -(shift+t+2499-SpikeLen/2)+1 ))=1;
    a=gausswin(SpikeLen);
    SpikeRate=conv(gausswin(SpikeLen),SpikeVec)/sum(a);
    SpikeRate(1:SpikeLen)=[];
    SpikeRate(ConvLen-5000+SpikeLen-SpikeLen+1:end)=[];
    
    block(1,X_ind')= block(1,X_ind')+SpikeRate;
    block(2,X_ind')= block(2,X_ind')+1;
end

removeFirst=find(sum(block)>0,1,'first');
removeLast=find(sum(block)>0,1,'last');

Y=block(1,:)./(block(2,:)+.000001);
X=[1:size(Y,2)]/10^2-6;

remove=find(block(2,:)<5);%Here take off noise 
Y(remove)=[];
X(remove)=[];

if cell==2 || cell==3
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

X_ind=[];FR=[];xSTA=[];nStim=[];
flag=1;ConvLen=10*10000;%500000;
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
r=1;row=[];
for lo=StimTime1(1:end-1)';
    First=find(TT(1,cell).sp>lo,1,'first');
    Last=find(TT(1,cell).sp<=lo+2000000,1,'last');
    AP(end+1:end+Last-First+1)=TT(1,cell).sp(First:Last)-lo;
    row(end+1:end+Last-First+1)=r;
    r=r+1;
end

nbin=20000;
SpikeLen=10;
figure(7)

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
figure

plot(eig(STC),'.')
xlim([0 200])
xlabel('Eigenvalue number');ylabel('Eigenvalue (Variance)');
title(['Cell num. ' num2str(cell)]);

[V,D]=eig(STC)

figure
plot(V(:,200));
title(['Cell num. ' num2str(cell)]);
