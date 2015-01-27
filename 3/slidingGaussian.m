function [y] = slidingGaussian(x, SLIDING_WINDOW, SIGMA)
    % does 1d gaussian filter on signal x, window size is SLIDING_WINDOW
    % returns a vector with the same size as x

        %window shall be pair, we divide it by 2
        if (mod(SLIDING_WINDOW,2)~=0)
           SLIDING_WINDOW = SLIDING_WINDOW+1;
        end
        w = linspace(-SLIDING_WINDOW / 2, SLIDING_WINDOW / 2, SLIDING_WINDOW);
        kernel = 1/(sqrt(2*pi)*SIGMA)*exp(-1.*(w .^ 2/(2*(SIGMA^2))));
        
        %{
        kernel2 = exp(-w .^ 2 / (2 * SIGMA ^ 2));
        kernel2 = kernel2 / sum (kernel2); % normalize
        %}

        
        y = conv(x, kernel, 'same');
        %y2 = conv(x, kernel2, 'same');
    end