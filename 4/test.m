%see http://en.wikipedia.org/wiki/Spike-triggered_covariance
%see Matlab Code in wiki: http://pillowlab.cps.utexas.edu/code_STC.html
%{
STC=0;
for i=1:length(accSTAStims)
    STC = STC + accAPs(i).*(accSTAStims(i,:)-STA')*(accSTAStims(i,:)-STA')';
end

STC = 1/(totalAPs-1).*STC

C=0;
for i=1:length(accSTAStims)
    C = C + accSTAStims(i,:)*accSTAStims(i,:)';
end

C=1/(totalStims-1).*C
%}

STC = accSTAStims'*(accSTAStims.*repmat(accAPs,1,STIMS_IN_STA_WINDOW)) ...
    / (totalAPs-1) - STA*STA'*totalAPs/(totalAPs-1);

[u,s,v] = svd(STC); 
%s is a diag matrix with the eigenvalues


%  Plot results -----------
figure;
plot(diag(s), 'o'); % examine eigenvalues
title('eigenvalues (using svd)');


[V,D] = eig(STC)
%produces matrices of eigenvalues (D, diag(D) are eigenvalues) 
%and eigenvectors (V, its columns are the eigenvectors of A) of matrix A, 
%so that A*V = V*D.

%  Plot results -----------
figure;
plot(diag(D), 'o'); % examine eigenvalues
title('eigenvalues (using eig)');
