clc;
%clear all;
close all;
%figure('units','normalized','outerposition',[0 0 1 1])
%{
figureEx('Custom', 'Maximize')

x = linspace(0,3*pi,200);
y = cos(x) + rand(1,200);
scatter(x,y)
%}

%{
actionPotentials = uint32(50:10:200);
IDX = uint32(1:1000);
res = zeros(length(actionPotentials), 3);
res(:,1) = actionPotentials;
res(:,2) = actionPotentials(:)>100 & actionPotentials(:)<150;
res(:,3) = IDX(1:length(actionPotentials));


ind = IDX(logical(res(:,2)))
%}


spikesAtTimes = zeros(75,1);
spikeTimes = [30;12;24;15;55;74;54;31;23;64;75];
%sortedSpikeTimes = sort(spikeTimes);
spikesAtTimes(spikeTimes)=1;

indexes = 1:75;
filter = spikesAtTimes(:,1) > 0;
spikeTimes2 = indexes(filter)';
accPerTime = spikesAtTimes(filter);

%accPerTime = [1;1 ;1 ;1 ;1;1;1;1;1;1;1;1];
%bins = 1:20:60;
bins = [0;10;25;50;75;90];


[bincounts,binIndex] = histc(spikeTimes,bins);
sumByBins = accumarray(binIndex,accPerTime, [size(bins,1) 1]);
figure;
h1 = bar(bins,sumByBins,'histc');
set(h1,'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0.60 0.60 0.60], ...
    'FaceAlpha',0.70);
hold on;
h2 = bar(bins,sumByBins/2,'histc');
set(h2,'FaceColor','k','EdgeColor','k','FaceAlpha',0.40);
    

%{
BIN_SIZE=100;
NUM_OF_BINS = STIMULUS_WINDOW/BIN_SIZE;
bins = 0:BIN_SIZE:STIMULUS_WINDOW;

indexesOfStimulusTime = 1:STIMULUS_WINDOW;
filter = accAPs(:,1) > 0;
spikeTimes = indexesOfStimulusTime(filter)';
accAPsAtLeastOne = accAPs(filter);




[bincounts,binIndex] = histc(spikeTimes,bins);
sumByBins = accumarray(binIndex,accAPsAtLeastOne, [length(bins) 1]);
figure;
bar(bins,sumByBins,'histc');
%}



%{ 
%light on/off
onOffData = zeros(STIMULUS_WINDOW,1);
onOffData(1:DELAY_TO_ON)=1;
h1=bar(onOffData,'y');
hold on;
onOffData(1:DELAY_TO_ON)=0;
onOffData(DELAY_TO_ON+1:STIMULUS_WINDOW)=1;
h2=bar(onOffData);
set(h2,'FaceColor','k','EdgeColor','k');

%get(bH,'Type') %For illustration purpose only.
pH = arrayfun(@(x) allchild(x),h1);
set(pH,'FaceAlpha',0.25);
%}
%{
set(h1,'FaceColor','r','EdgeColor','r');
hold on;
h2 = bar([1 2],[-0.0001 1.00001],'histc');
%set(h2,'FaceColor','y','EdgeColor','y');
%}
%{
% transparency smaple
axis([0 STIMULUS_WINDOW 0 1])

figure;
bH = bar(randn(10)); 
%get(bH,'Type') %For illustration purpose only.
pH = arrayfun(@(x) allchild(x),bH);
set(pH,'FaceAlpha',0.25);
%}