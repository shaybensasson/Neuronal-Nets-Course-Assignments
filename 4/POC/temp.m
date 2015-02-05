clc;
accBinned = accBinnedOrig;
%accBinnedOrig = accBinned;

accBinned_sorted = sortrows(accBinned, 1);

times = accBinned(:,1);
stimsAfterLinearFilter = accBinned(:,2);
aps = accBinned(:,3);

bins = 0:curBinSize:TICKS_IN_WINDOW;
[~,binIndex] = histc(times(:,1),bins);
%insert any unbinned data to last bin - WE SHALL HAVE NO UNBINNED
%binIndex(binIndex==0)=max(binIndex);

%group by mean
grp_aps = accumarray(binIndex, aps, [length(bins) 1], @mean, NaN);
grp_stimsAfterLinearFilter = accumarray(binIndex, stimsAfterLinearFilter, [length(bins) 1], @mean, NaN);

rvsrest = [ceil(bins') grp_aps grp_stimsAfterLinearFilter NaN(length(bins), 1)];

rvsrest(isnan(rvsrest(:,2)) | isnan(rvsrest(:,3)), :) = [];
grp_stimsAfterGenerator = curNeuron.PSTH{iBinSize}.Generator(rvsrest(:,3));
rvsrest(:,4) = normalize(grp_stimsAfterGenerator, NORMALIZE, NORMALIZE_BY_MAX);




%{
BIN_SIZE = 2
windowData = [1,2;2,4;3,6;4,8;5,10;6,12;7,14;8,2;9,1;10,NaN;10,5]
%tempWinData = [1,1;2,2;3,3;3,9];

times = windowData(:, 1);
vals = windowData(:, 2);
bins = times(1):BIN_SIZE:times(end);
[bincounts,binIndex] = histc(times(:,1),bins);
%insert any unbinned data to last bin
binIndex(binIndex==0)=max(binIndex);

%group by mean
%meanIgnoreNaNs = @(vector) mean(vector(~isnan(vector(:))));
sumIgnoreNaNs = @(vector) sum(vector(~isnan(vector(:))));
grouped = accumarray(binIndex, vals, [length(bins) 1], sumIgnoreNaNs)
%}