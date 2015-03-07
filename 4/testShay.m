TICKS_PER_BIN = 100;
nbins=REP_TRAIL_IN_AP_TICKS/TICKS_PER_BIN; %100 ticks each bin
SpikeLen=2; %how many spikes we consider 1 spike (smoothing)
figure;

[binCounts,xCenters]=hist(APTimesCrossTrails/AP_TICKS_PER_SECOND,nbins);
hold on;

gaussWinFilter=gausswin(SpikeLen);
yGauss=conv(gaussWinFilter,binCounts/countTrails)/sum(gaussWinFilter);%units:1/sample
%lowpassfilter: start after half a filter window and end before half
plot(xCenters,binCounts/countTrails); %TODO: remove
plot(xCenters,yGauss(SpikeLen/2:end-SpikeLen/2),'r'); 

ylabel('FR [10,000*Hz]');
xlim([t/AP_TICKS_PER_SECOND t/AP_TICKS_PER_SECOND+10]);