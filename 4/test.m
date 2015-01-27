clear all;
clc;

%% filter range of values
rawData = 1900:5:2100;

from = 1982;
to_excluded = 2015;

indexes = 1:length(rawData);
%declare indexes for the whole data range

filter = logical(rawData(:) >= from & rawData(:) < to_excluded);
%contains 1 where we care

indexesOfFilteredData = indexes(filter);
%get all the indexes that have 1s

filteredData = rawData(indexesOfFilteredData);
%filter the data by these indexes