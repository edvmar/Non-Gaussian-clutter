function [pFA, pTD] = NeymannPearsonLRT(r, toeplitzMatrixInverse, signal, treshold, h_N)
    
    q0 = r'*toeplitzMatrixInverse*r;
    q1 = (r-signal)'*toeplitzMatrixInverse*(r-signal);
    
    if h_N(q1)/h_N(q0) > treshold
        
    end
    
end

