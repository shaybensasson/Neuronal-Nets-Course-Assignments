xlim([1,5]);
hold('on');
h1a = plot(1:5,     11:15, '.-', 'LineWidth',1, 'DisplayName',' 0.5');
h1b = plot(1.5:5.5, 11:15, '.-', 'LineWidth',1, 'DisplayName',' 1.0', 'Color',h1a.Color);  % 100% opaque
h1a.Color(4) = 0.5;  % 50% transparent
h2a = plot(3:7,  15:-1:11, '.-r', 'LineWidth',1, 'DisplayName',' 0.3'); h2a.Color(4)=0.3;  % 70% transparent
h2b = plot(2:6,  15:-1:11, '.-r', 'LineWidth',1, 'DisplayName',' 0.7'); h2b.Color(4)=0.7;  % 30% transparent
h2c = plot(1:5,  15:-1:11, '.-r', 'LineWidth',1, 'DisplayName',' 1.0');  % 100% opaque = 0% transparent
legend('show','Location','west')