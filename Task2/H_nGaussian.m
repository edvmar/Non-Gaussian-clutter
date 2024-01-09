%%%%%%%%%%%%%%%%%%%% H_nGaussian %%%%%%%%%%%%%%%%%%%%%%%
%
% Tail distribution for the complex Gaussian distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h_N = H_nGaussian(y, numberOfPulses, sigma)
    
    N = numberOfPulses;
    eta = sigma^2;
    h_N = exp(-y/eta) / eta^N;

end

