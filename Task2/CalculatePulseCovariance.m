%%%%%%%%%%% CalculatePulseCovariance %%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates Sigma as a Toeplitz matrix and L : Sigma = LL^T
% by Cholesky decomposition with diagonal loading.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = CalculatePulseCovariance(numberOfPulses, epsilon, delta)
    
    sigma = eye(numberOfPulses);
    
    for i = 1:numberOfPulses
        for j = 1:numberOfPulses
            sigma(i,j) = exp(-(i-j)^2*delta);
        end
    end
    
    L = chol(sigma + epsilon*eye(numberOfPulses));
end

