function [ y ] = funLinearWindowWithZeros( x, h )
%FUNLINEARWINDOW passes a linear window h on the input x
%   Expects x and h to be row vectors
    %see http://matlabtricks.com/post-3/the-basics-of-convolution
    %y = zeros(length(h)-1 + length(x), 1);
    
    xPad = x;
    xPad = padarray(xPad, [length(h)-1 0], 0, 'pre');
    xPad = padarray(xPad, [length(h)-1 0], 0, 'post');
    
    yPad = zeros(length(h)-1 + length(x), 1);
        
    for i = 1:(length(h)-1 + length(x))
        xx = xPad(i:i+length(h)-1);
        yPad(i, 1) = h'*xx;
    end
    
    yPad(length(yPad)-(length(h)-1)+1:length(yPad)) = [];
    yPad(1:(length(h)-1)) = [];
    
    y = yPad;
end

