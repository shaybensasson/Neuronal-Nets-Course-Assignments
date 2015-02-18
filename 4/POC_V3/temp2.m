figure;
f1 = STA;
f1 = normalize(f1, 1, 1);

f2 = Filter;
f2 = normalize(f2, 1, 1);

x = linspace(-Simulation.STA_WINDOW_IN_MS, 0, length(f1));
plot(x, f1);

hold on

plot(x, f2);

%ylim([minSTAValue-0.1 maxSTAValue+0.1]);


title(sprintf('Neuron #%d', iNeuron));
xlabel('Time (ms)');
ylabel('Light levels');

CreateTitleForSubplots('\bf STA');