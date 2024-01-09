%%%%%%%%%%%% Calculate q1 with ML estimate of alpha %%%%%%%%%%%%%
%
% Calculates q1 by eq (7.34) in Modern Radar Detection by Greco/de Maio
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function q1 = CalculateQ1WithEstimatedAlpha(z,toeplitzMatrixInverse,...
                                            steeringVector, steeringVectorNorm)

    p = steeringVector;

    term1 = MultidimensionalNorm(z, toeplitzMatrixInverse);

    q1 = term1 - abs(p'*toeplitzMatrixInverse*z).^2/steeringVectorNorm;
    
    q1 = real(q1);

        
end

