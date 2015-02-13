%find times>=40*Simulation.TICKS_IN_SECOND
onset = times(1,1)+40*Simulation.TICKS_IN_SECOND;
[row,~] = find(times(:,1) >= onset);
I2 = times(times(:,1)>=onset);
I2(1,1)