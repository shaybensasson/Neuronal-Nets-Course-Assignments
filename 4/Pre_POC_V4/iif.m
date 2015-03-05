function out = iif(cond,a,b)
%IIF implements a ternary operator

% pre-assign out
if (cond)
    out = a;
else
    out = b;
end