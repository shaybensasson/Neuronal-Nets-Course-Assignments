%{
x = [0:.1:200];
norm = normpdf(x,100,7);
%}
%{
x = [-100:.1:100];
norm = normpdf(x,0,1);

figure;
plot(x,norm);
%}

Mu = 10
Sigma = 30
I = normrnd(Mu,Sigma)