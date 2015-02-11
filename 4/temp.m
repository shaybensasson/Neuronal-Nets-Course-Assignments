x=[0:1:10];
figure
subplot(1,2,1), plot(x,x.^2)
subplot(1,2,2), plot(x,x.^3-1)


% and then, type:

%ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text  = '\bf Do you like this title?';
CreateTitleForSubplots(text);