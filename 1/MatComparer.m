clear all;
load('single.outArrSingle.mat');
load('multi.outArr.mat');


figure(1);
hold on;
comp = 1;
VmComp1 = outArr(:,comp,1);
plot(VmComp1, '-b');
comp = 2;
%Vm
VmComp2 = outArr(:,comp,1);
plot(VmComp2, '--g');
hold off;

figure(2);
hold on;
%n
plot(outArr(:,comp,2), '--r');
%m
plot(outArr(:,comp,3), '--g');
%h
plot(outArr(:,comp,4), '--b');

hold off;

%% cut by Param
figure(3);
hold on
VmForAllComps = outArr(:,:,1); %seperate Vm from other params
comp=1;
VmComp1=VmForAllComps(:,comp);
plot(VmComp1, '-k');

plot(VmComp2, '--y'); %from upstairs
hold off;