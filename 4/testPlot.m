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
xFrom = t/AP_TICKS_PER_SECOND;
SECONDS_TO_DISPLAY = 10;
xTo = t/AP_TICKS_PER_SECOND+SECONDS_TO_DISPLAY;
xlim([xFrom xTo]);

hold on
shift=length(Filter)/2; %because it's a low pass filter, we must shift by half and reduce half at the end
estRate = SimulatedFR*100;
plot((0+(t+shift:t+TRAIL_IN_SIM_TICKS-(shift+1)))/AP_TICKS_PER_SECOND,estRate,'k');
xlabel('Time [sec]');
title(['Cell num. ' num2str(iNeuron)]);
legend('PSTH','Simulated FR')

valFrom = t/AP_TICKS_PER_SECOND;
actualRate = actualRate(xCenters>=xFrom & xCenters<=xTo);
rates = [actualRate estRate];
actualRate = rates(:,1);
estRate = rates(:,2);
Diff = actualRate-estRate;
SSE = norm(Diff,2)^2; %lower is less err
Cor = abs(xcorr(actualRate,estRate,0,'coeff')); % 1 if are equal

PositionOfSubplot = get(gca, 'Position');

ha = axes('Position',PositionOfSubplot,'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        
text(0, 0.05, ...
    sprintf(['Diff: %.2f;\n' ...
    'SSE: %.2f;\n' ...
    'Cor: %.4f'], ...
        Diff, SSE, Cor), ...
'HorizontalAlignment' ,'left','VerticalAlignment', 'bottom');
