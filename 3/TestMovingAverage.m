function TestMovingAverage() 
    %see: http://matlabtricks.com/post-11/moving-average-by-convolution

    x = [1 7 1 4 4 7 1];
    SLIDING_WINDOW = 3;
    format rat
    movingAverage(x, SLIDING_WINDOW)
end