function correctPrecent = GetCorrectPrecentage(matrix1, matrix2)
    matrixcomp=sum(matrix1==matrix2);
    totalNodes = size(matrix1,1) *size(matrix1,2);
    if (matrixcomp~=0)
    correctPrecent = matrixcomp/totalNodes;
    else
        correctPrecent=0;
    end
    
end
