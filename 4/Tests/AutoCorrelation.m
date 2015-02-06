% Authors
%   Reza Arfa <rezaarfa (at)gmail.com>
%
% License
%   The program is free to use for non-commercial and academic purposes.
%
% Changes
%   01/01/2013 (!)

clc; close all; clear all;

% ========== init some variables:
iteration = 200;
fs = 10000;
t = -1:1/fs:1;

iterationSteps   = 1:1:iteration;
autoCorrelation = zeros(1,iteration);
convolution      = zeros(1,iteration);

set(gcf,'Color',[1 1 1]);


% ========== plot the same signal twice (y1 and y2):
disp('Let y1=y2 to be defined as following figures:')
disp('============================================================')
w=1;
%f = @(t,w) rectpuls(t,w);
f = @(t,w) tripuls(t,0.5,-1);
y1 = f(t,w);
y2 = f(t,w);
subplot(4,2,1);plot(t,y1,'Color','blue','LineWidth',2),axis([-1 1 -0.2 1.2]);
ylabel('y1')
subplot(4,2,2);plot(t,y2,'Color','red','LineWidth',2),axis([-1 1 -0.2 1.2]);
ylabel('y2')


% ========== Auto Correlation of a signal (y1=y2):
disp(' ');disp(' ')
disp('Auto Correlation of a signal (y1=y2):')
disp('============================================================')

for i = 1:iteration
    
    moveStep = (i-100)/100;
    y1 = f(t,1);
    y2 = f(t-moveStep,1);
    
    subplot(4,2,3:4)
    hold off;plot(t,y1,'Color','blue','LineWidth',2),axis([-1 1 -0.2 1.2]);
    hold on; plot(t,y2,'Color','red', 'LineWidth',2),axis([-1 1 -0.2 1.2]);
    
    
    autoCorrelation(i) = sum(y1.*y2);
    
    subplot(4,2,5:6)
    hold off
    plot(iterationSteps(1:i),autoCorrelation(1:i),'Color','black','LineWidth',2); 
    axis([1 iteration -100 fs+100]);
    xlabel('t')
    ylabel('AutoCorrelation(y1, y2)(t) ')
    
    pause(0.01)
end
