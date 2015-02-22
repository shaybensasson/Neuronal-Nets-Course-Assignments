x = (0:100)';
y = 5 + 10./(1+exp(-(x-40)/10)) + randn(size(x));
plot(x,y,'bo')

%%
f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));
p = nlinfit(x,y,f,[0 20 50 5])
line(x,f(p,x),'color','r')

%%
fit(x,y,'a + b ./ (1 + exp(-(x-m)/s))','start',[0 20 50 5])
%%
%idxs = 1:length(m);
%m = [m idxs'];

filter = m(:,1)>=-0.6 & m(:,1)<=0.6;
mFiltered = m(filter, :)
idxFrom = mFiltered(1,4);
m(1:idxFrom, 2) = m(idxFrom, 2);
idxTo = mFiltered(end,4);
m(idxTo:end, 2) = m(idxTo, 2);

xData = m(:, 1);
yData = m(:, 2);