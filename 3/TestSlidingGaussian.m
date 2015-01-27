function TestSlidingGaussian() 
    %x = [1 1 1 1 1 1 1];
    %x = rand(100,1);
    x = 1:100;
    SLIDING_WINDOW = 5;
    SIGMA = 1;
    [y] = slidingGaussian(x, SLIDING_WINDOW, SIGMA)
    
    figure;
    plot(x, 'r');
    hold on;
    plot(y, 'b');
    
    %plot(y2, 'g');
    %axis([0 100 0 1]);
    
end