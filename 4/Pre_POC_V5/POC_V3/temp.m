x = funcXData;
y1 = funcYMeans;
y2 = bincounts;

%{
x1 = 0:0.1:40;
y1 = 4.*cos(x1)./(x1+2);
x2 = 1:0.2:20;
y2 = x2.^2./x2.^3;
%}

figure
%line(x,y1,'Color','r')
h = plot( fitresult, x, y1 );

%{

ax1 = gca; % current axes

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none');
ax2.YColor = [100 100 100] / 150;
ylabel(ax2, 'Bin counts');

hl = line(x,y2,'Parent',ax2,'Color','k');
hl.LineStyle = '--';
hl.Color(4) = 0.5;

grid on;
%}

ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none');
ax2.XColor = [0 0 0 0]; %transparent
ax2.YColor = [100 100 100] / 150;
ylabel(ax2, 'Bin counts');

hl = line(x,y2,'Parent',ax2,'Color','k');
hl.LineStyle = '--';
hl.Color(4) = 0.5;

grid on
