function [ data ] = remap( data, OldMin, OldMax, NewMin, NewMax)
%REMAP Summary of this function goes here
%   Detailed explanation goes here
    OldRange = (OldMax - OldMin);
    NewRange = (NewMax - NewMin);
    data = (((data - OldMin) * NewRange) / OldRange) + NewMin;

end

