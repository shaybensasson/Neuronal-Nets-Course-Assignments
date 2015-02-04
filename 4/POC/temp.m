clc;
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
