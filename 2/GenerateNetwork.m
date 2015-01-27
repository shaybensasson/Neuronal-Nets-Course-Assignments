function net = GenerateNetwork(n)

A=round(rand(n,n));
A(A==0)=-1;
net = A;


end