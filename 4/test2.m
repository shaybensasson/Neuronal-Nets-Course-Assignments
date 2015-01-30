%A=[5,8,16;4,1,8;-4,-4,-11]; %real eigvects

%http://math.mit.edu/~gs/linearalgebra/ila0601.pdf
A=[.8,.3;.2,.7];
[V,D] = eig(A);
[evals,idx] = sort(diag(D),'descend');
EigenValues = evals
evects = V(:,idx); 
EigenVects = evects

%{
EigenValues =

          0.5
            1


EigenVects =

     -0.70711      0.83205
      0.70711       0.5547
%}