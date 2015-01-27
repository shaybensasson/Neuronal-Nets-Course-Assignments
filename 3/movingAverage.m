function y = movingAverage(x, SLIDING_WINDOW)
% does moving average on signal x, window size is SLIDING_WINDOW
% returns a vector with the same size as x


   %MA assumes the window is odd, for symmetry
   if (mod(SLIDING_WINDOW,2)==0)
       SLIDING_WINDOW = SLIDING_WINDOW+1;
   end
   kernel = ones(1, SLIDING_WINDOW) / SLIDING_WINDOW;
   y = conv(x, kernel, 'same');
end