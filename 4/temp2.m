figure('Name', 'STA');
plot(normalize(STA, 1, 1));

figure('Name', 'STAPerTrail');
for idx=1:10
%for idx=1:size(STAPerTrail,1)
    hold on;
    tData = normalize(STAPerTrail(idx, :), 1, 1);
    plot(tData);
end