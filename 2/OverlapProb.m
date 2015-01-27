%BUILD A HOPFIELD NETWORK WHICH RETRIEVES RANDOM MEMORIES
clear all; %clear all variables
close all; %close all open figure windows

%DEFINE NETWORK PARAMETERS
NumPatt = 15 %Number of memory patterns that network will attempt to store and retrieve
N = 50 %Number of neurons in the network (and number of elements in a memory pattern)
NumTimeSteps = 50 %length of run

%MEMORY PATTERNS
Mem_mat = 2*round(rand(N,NumPatt))-1; %random strings of 1's and -1's
%{
Note that the round command rounds off the random number, giving a 1 with probability 0.5
(i.e. if 0.5 <= rand < 1) or a 0 with probability 0.5 (if 0 < rand < 0.5). Multiplying 0 or 1 by 2
and then subtracting 1 gives either -1 or +1, respectively. Thus, this line defines a matrix whose
NumPatt columns each contain an N-element vector of random memory patterns. 
%}

%DEFINE SYNAPTIC WEIGHT MATRIX - Hebb Rule used to define synaptic weights
W_mat = Mem_mat*Mem_mat';

%INITIALIZE NETWORK ACTIVITY
%create an input pattern from the 1st pattern with noise
x_vect=zeros(1,N);
qStart = 0.1; %average overlap of initial activity pattern with first memory (the one to be retrieved)
selectedPattern = randi(NumPatt, 1,1)
for i=1:N
 if rand < qStart %assign x_i equal to Mem(i,1)
 x_vect(i)=Mem_mat(i,selectedPattern);
 else %assign random value for x_i
 if rand < .5
 x_vect(i) = 1;
 else
 x_vect(i) = -1;
 end
 end
end

initial_x_vect = x_vect;
q_mat = zeros(NumPatt, NumTimeSteps);
for i=1:NumPatt
    q_mat(i,1)= x_vect*Mem_mat(:,i)/N;
end


figure(1)
num_of_plots = (NumTimeSteps/10)+1;
num_of_plots = num_of_plots *2; %x and q at step

tStep=1;
iPlot=1;
h(iPlot) = subplot(num_of_plots,1,iPlot);
image(50*x_vect)
xlabel(strcat('x(t=',num2str(0),')'))
ylabel('blue=-1,red=+1')

iPlot=iPlot+1;
h(iPlot) = subplot(num_of_plots,1,iPlot);
%plot(abs(q_mat(:,tStep)'), '.');
%xlim([0 NumPatt+1]);
%ylim([0 1]);
colormap('hot')
imagesc(abs(q_mat(:,tStep)'))
colorbar
xlabel(strcat('q(t=',num2str(0),')'))


%RUN MODEL AND Show X Vector at each step

for tStep=2:NumTimeSteps
    x_vect = (2*(W_mat*x_vect' >= 0) - 1)'; %update rule
    for i=1:NumPatt
       q_mat(i, tStep)= x_vect*Mem_mat(:,i)/N;
    end
    
    if (mod(tStep,10)==0)
         iPlot=iPlot+1;

         h(iPlot) = subplot(num_of_plots,1,iPlot);
         image(50*x_vect)
         xlabel(strcat('x(t=',num2str(tStep-1),')'))

         iPlot=iPlot+1;
         h(iPlot) = subplot(num_of_plots,1,iPlot);
         %plot(abs(q_mat(:,tStep)'), '.');
         %xlim([0 NumPatt+1]);
         %ylim([0 1]);
         
         colormap('hot')
         imagesc(abs(q_mat(:,tStep)'))
         colorbar
         
         xlabel(strcat('q(t=',num2str(tStep-1),')'))

         %ylabel('blue=-1,red=+1')
     end
end



figure(2)
subplot(3,1,1)
image(50*initial_x_vect)
%ylabel('initial x')
xlabel(strcat('x(t=',num2str(0),')'))

subplot(3,1,2)
image(50*x_vect)
%ylabel('final x')
xlabel(strcat('x(t=',num2str(NumTimeSteps-1),')'))

subplot(3,1,3);
%{
dif = Mem_mat(:,1)'-x_vect;
image(50*dif)
xlabel(strcat(num2str(length(find(dif>0))), ' bits were flipped'))
ylabel('red=diff')
%}

diff_vs_all_others = Mem_mat-repmat(x_vect',1,NumPatt);

image(50*diff_vs_all_others)
xlabel('stored patterns');
ylabel('neuron')