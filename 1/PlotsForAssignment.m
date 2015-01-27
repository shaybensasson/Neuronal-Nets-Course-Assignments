clc;

%% Q5

x = Acc(:,1);
y = Acc(:,2);
a = 50;

cmp = ones(size(x)); % create the color maps changed as in jet color map
scatter(x, y, 10, cmp, 'filled');

xlabel('External Current');
ylabel('# of Action Potentials');
title('Q6');

