clc;
clear all;
tic
pause
t=toc;
disp(['Elapsed time is ' datestr(datenum(0,0,0,0,0,t),'HH:MM:SS') '.']);
return;
%{
N = 32;
im = [1 2;3 4];
W = kron(im,im);
W
%{
% Erasing self weight
ind = 1:N^2;
f = find(mod(ind,N+1)==1);
W(ind(f),ind(f)) = 0;
%}

W2 = GetWeights(im)
%}


net = magic(5)
net = circshift(net, [0, 2])
return;
folder = '.\images\cd01A\';
dir([folder '*.jpg'])
return;

noiseNumber = 5;
NodesPerPattern = 10;

clc
%for i=1:5
    clc
    %p = randperm(5)
    %ind= ceil(rand(noiseNumber,1)*(NodesPerPattern-1)));
    p = randperm(NodesPerPattern)
    p(1:noiseNumber)
    %{
    ind= round(rand(noiseNumber,1)*(NodesPerPattern-1))+1
    is = sum(mod(ind,1))>0;
    %}
%end

%{
curNoise=0:0.05:1;
ceil(NodesPerPattern*curNoise)
%}


%mod(1.5,1)
%mod(40,1)