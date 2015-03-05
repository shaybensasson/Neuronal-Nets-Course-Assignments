map = [43, 87, 154; ... %blue
    32, 162, 58; ... %green
    0, 0, 0; ... % black
    207, 210, 146; ... %yello for stims
    205, 151, 151]; %red for aps

map = map./255;

h = plot(1:10,1:10);
h.Color = map(5,:);