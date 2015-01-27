function result = iif(condition,trueExprStr,falseExprStr)
   narginchk(3,3);
   if condition
     result = trueExprStr;
   else
     result = falseExprStr;
   end
end