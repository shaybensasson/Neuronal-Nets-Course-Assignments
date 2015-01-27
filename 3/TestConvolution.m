%see: http://matlabtricks.com/post-3/the-basics-of-convolution

clc;


x = [1, 2, 1, 3];
h = [2, 0, 1];
y = conv(x, h)

%Convolution is commutative, the order of the operands makes no difference:
yComm = conv(h, x)

%same length as input (supressing edges)
ySame = conv(x, h, 'same') 

%all the conv result
yFull = conv(x, h, 'full')

%Option valid returns those elements only which were fully covered, 
% so there was no sliding off during the windowing
yValid = conv(x, h, 'valid')
