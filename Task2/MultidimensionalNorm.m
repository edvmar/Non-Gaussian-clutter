function xprimMx = MultidimensionalNorm(vectorMatrix,matrix)
    matrix2 = matrix*vectorMatrix;
    xprimMx = dot(vectorMatrix,matrix2);
end

