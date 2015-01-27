clc;

%% NvsP
%{
xspan=1:1:length(numOfCorrectInferences);
yspan = numOfCorrectInferences;

plot(xspan,yspan, 'b');
xlabel('Patterns Stored');
ylabel('Patterns Infered');
title(['NotNaturalImages' ' - Stored vs Infered' ' - ' ...
    iif(SYNC_UPDATE, '', 'A') ...
    'Sync Update']);
%}


%% NvsP Iterations
%totalPatterns = length(numOfCorrectInferences);
totalPatterns = 70;
xspan=1:totalPatterns;

[ax,p1,p2] = plotyy(xspan,numOfCorrectInferences(1:totalPatterns), ...
    xspan, numOfIterationsUntilInf(1:totalPatterns), ...
    'plot', 'semilogy');

ylabel(ax(1),'Correct inferences', 'fontweight','bold','fontsize',14); % label left y-axis
ylabel(ax(2),'Total iterations afterall patterns were tested', ...
    'fontweight','bold','fontsize',12); % label right y-axis

set(ax(1), 'YLim', [1 totalPatterns], 'YTick', [1 totalPatterns]);
set(ax(2), 'YLim', [1e0 1e4], 'YTick', [1e0 1e1 1e2 1e3 1e4]);
xlabel('Patterns Stored', 'fontweight','bold','fontsize',14);




title(['NotNaturalImages' ' - Correct Inferences and Total Iterations' ' - ' ...
    iif(SYNC_UPDATE, '', 'A') ...
    'Sync Update'], ...
    'fontweight','bold'); 


%% Det/Random Noise
%{
numOfCorrectInferences(numOfCorrectInferences==0) = nan;
j = flipud(jet);
j(1,:) = [ 1 1 1 ];
colormap(j);
numOfPatterns=zeros(size(numOfCorrectInferences));
for iRow=1:size(numOfPatterns,1)
    numOfPatterns(iRow,:)=1:1:length(numOfPatterns);
end

rel = numOfCorrectInferences./numOfPatterns;

imagesc(rel);
colorbar;

xlabel('Patterns Stored');
ylabel('Noise');
title(['NotNaturalImages' ' - ' ...
    iif(RANDOM_NOISE, 'Random', 'Deterministic') ...
    ' Noise ' '- ' ...
    iif(SYNC_UPDATE, '', 'A') ...
    'Sync Update']);

set(gca,'YTick',1:3)
set(gca,'YTickLabel', ['30%';'40%';'50%'])

%imagesc(numOfCorrectInferences);
%}

%% Natural Images/Diverted
%{
xspan=1:1:length(numOfCorrectInferences);
yspan = numOfCorrectInferences;

plot(xspan,yspan, 'b');
xlabel('Patterns Stored');
xlim([0 length(numOfCorrectInferences)]);
ylim([0 length(numOfCorrectInferences)]);
ylabel('Patterns Infered');
title(['NaturalImages' ...
    iif(DIVERT, ' - Shifted - ', ' - ') ...
    'Stored vs Infered - ' ...
    iif(SYNC_UPDATE, '', 'A') ...
    'Sync Update']);
%}