y = [1 2 3];
yhat = [2 1 1];

%see http://www.mathworks.com/matlabcentral/answers/104189-calculate-sum-of-square-error
SSE = norm(y-yhat,2)^2; %lower is less err

%the similarity between two signals, we only need zero lag
CC = xcorr(y,yhat,0,'coeff'); % 1 if are equal