%%%%%%%%%%% CalculatePulseCovariance %%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the covariance as a Toeplitz matrix with unit diagonal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta)
    
    % Matrices for row and column indeces
    indeces = 1:numberOfPulses;
    iRow    = repmat(indeces',1, width(indeces));
    jColumn = repmat(indeces, width(indeces),1);
   
    toeplitzMatrix = exp(-(iRow-jColumn).^2*delta);
    
end

