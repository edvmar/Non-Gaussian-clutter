function h_N = H_nGaussian(y, numberOfPulses, sigma)
    
    % Tail distribution with derivatives
    N = numberOfPulses;
    eta = sigma^2;
    h_N = exp(-y/eta) / eta^N;

end

