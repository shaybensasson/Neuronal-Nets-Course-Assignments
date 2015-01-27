function W= GetWeights(net)

[p,q] = meshgrid(net(:), net(:));
pairs = p(:).*q(:);

n=size(net(:),1);

W = reshape(pairs,n,n);

%disp(['isequal? ' num2str(isequal(net'*net, W))]);

end