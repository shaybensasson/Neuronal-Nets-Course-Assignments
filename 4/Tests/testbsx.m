fun = @(A,B) nansum([A B], 2);
A = [NaN;2]
B = [10;20]

C = bsxfun(fun,A,B)