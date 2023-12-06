%%%%%%%%%%% CalculatePulseCovariance %%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates Sigma as a Toeplitz matrix and L : Sigma = LL^T
% by Cholesky decomposition with diagonal loading.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta)
    
    % Matrices for row and column indeces
    indeces = 1:numberOfPulses;
    iRow    = repmat(indeces',1, width(indeces));
    jColumn = repmat(indeces, width(indeces),1);
    
    % Toeplitz matrix
    toeplitzMatrix = exp(-(iRow-jColumn).^2*delta);
    
end

