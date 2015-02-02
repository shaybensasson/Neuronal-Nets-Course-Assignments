clc;
close all;

%{ 
linear filter the stimVals by the STA filter
we flip the STA to fix the convolution direction, hence,
moving the STA window from the start of the source to the end

see http://en.wikipedia.org/wiki/Cross-correlation
%}
%stimsAfterLinearFilter = conv(stimValuesOnWindow, flipud(STA'), 'same');

%source = 10:10:50
%source = ones(10,1);
source = round(sin(1:200)*10);
%filter = ones(3,1);
source = source';

%filter = ones(2,1);
%filter = [1 2 -1 2]';
filter = STA;
%filter = [1 2 2 1]';

resultValid = conv(source, flipud(filter'), 'valid')
resultSame = conv(source, flipud(filter'), 'same')
resultFull = conv(source, flipud(filter'), 'full')

figure;
subplot(3,1,1);
y = source;
x = 1:length(y);
plot(x,y,'k');
strValues = strtrim(cellstr(num2str(y,'(%d)')));
text(x, y,strValues,'VerticalAlignment','top');

subplot(3,1,2);
y = filter;
x = 1:length(y);
plot(x,y,'b');
strValues = strtrim(cellstr(num2str(y,'(%.2f)')));
text(x, y,strValues,'VerticalAlignment','top');
%set(gca,'XTick',x)

subplot(3,1,3);
y = resultValid;
%y(end-length(filter)+1+1:end, 1) = NaN; %discard non full filtered points
x = 1:length(y);
plot(x,y,'r');
strValues = strtrim(cellstr(num2str(y,'(%.2f)')));
text(x, y,strValues,'VerticalAlignment','top');
%set(gca,'XTick',x)
