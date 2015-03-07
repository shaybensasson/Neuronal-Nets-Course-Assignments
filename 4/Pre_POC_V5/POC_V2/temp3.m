bins = 1:Simulation.STIMULUS_EACH_TICKS:STA_WINDOW_IN_TICKS;
[bincounts,binIndex] = histc(1:STA_WINDOW_IN_TICKS,bins);

%insert any unbinned data to last bin
unbinnedIndex = max(binIndex);
%last bin was not created
if (~bincounts(end))
    unbinnedIndex = unbinnedIndex+1;
end
binIndex(binIndex==0)=unbinnedIndex;

%group by times and mean stim vals
meanIgnoreNaNs = @(vector) mean(vector(~isnan(vector(:))));
grp_mean = accumarray(binIndex', accSTAWindow', [length(bins) 1], meanIgnoreNaNs, NaN);

binned = [bins' grp_mean];

plot(binned(:,1), binned(:,2));
    