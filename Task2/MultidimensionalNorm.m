%%%%%%%%%%% MultidimensionalNorm %%%%%%%%%%%%%%%%%%%%%%%
%
% Performs the calculation of the quadratrue x'Ax
% Where x as an extra dimension of sample sizes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xprimMx = MultidimensionalNorm(vectorMatrix, matrix)

    matrix2 = matrix*vectorMatrix;
    xprimMx = dot(vectorMatrix,matrix2);

end

