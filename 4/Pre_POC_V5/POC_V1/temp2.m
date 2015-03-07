%http://www.mathworks.com/help/signal/ug/residual-analysis-with-autocorrelation.html
data = curNeuron.PSTH{iBinSize}.RVsRest;
times = data(:,1);
PSTH = data(:,2);
        
stimsAfterLinearFilter = data(:,3);
stimsAfterGenerator = data(:,4);
        
%x = -3:0.01:3;
x = times;
%rng default
%y = 2*x+randn(size(x));
y = PSTH;
plot(x,y)


coeffs = polyfit(x,y,1);
yfit = coeffs(2)+coeffs(1)*x;

plot(x,y)
hold on
plot(x,yfit,'linewidth',2)


residuals = y - yfit;
%[xc,lags] = xcorr(residuals,50,'coeff');
lag = ceil(TICKS_IN_WINDOW/PSTH_BIN_SIZES(iBinSize));
[xc,lags] = xcorr(residuals,lag,'coeff');


CONF_INTERVAL = .95;
conf = sqrt(2)*erfcinv(2*(1-CONF_INTERVAL)/2);
%lconf = -conf/sqrt(length(x));
lconf = -(1-CONF_INTERVAL);
%upconf = conf/sqrt(length(x));
upconf = (1-CONF_INTERVAL);


figure
stem(lags,xc,'filled')
ylim([lconf-0.03 1.05])
hold on
plot(lags,lconf*ones(size(lags)),'r','linewidth',2)
plot(lags,upconf*ones(size(lags)),'r','linewidth',2)
title(sprintf('Sample Autocorrelation with %d%s Confidence Intervals', CONF_INTERVAL*100, '%'));

