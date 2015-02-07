function [ result ] = normalize(a, NORMALIZE, NORMALIZE_BY_MAX);
%NORMALIZE Summary of this function goes here
%   Detailed explanation goes here
    result = a;
    
    filter = ~isnan(a(:));
    if (NORMALIZE) 
        m = mean(a(filter));
        
        result = a-m;
                
        %normalize
        if (NORMALIZE_BY_MAX) 
            result = result/max(abs(result(filter)));
        end
    end
end

