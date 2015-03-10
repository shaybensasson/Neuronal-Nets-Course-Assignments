function [ y ] = funLinearWindow( x, h )
%FUNLINEARWINDOW passes a linear window h on the input x
%   Expects x and h to be row vectors
%   Nans must be excluded from all inputs
    %see http://matlabtricks.com/post-3/the-basics-of-convolution
    
    %we use NaNs to filter partial convolved values
    xPad = x;
    xPad = padarray(xPad, [length(h)-1 0], NaN, 'pre');
    xPad = padarray(xPad, [length(h)-1 0], NaN, 'post');
    
    yPad = zeros(length(h)-1 + length(x), 1);
        
    for i = 1:(length(h)-1 + length(x))
        xx = xPad(i:i+length(h)-1);
        yPad(i, 1) = h'*xx;
    end
    
    yPad(isnan(yPad)) = [];
    
    y = yPad;
end

