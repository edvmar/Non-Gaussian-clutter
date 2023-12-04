function h_N = TailDistributionComplexGaussian(x, N, sigma)
    
    % Tail distribution with derivatives
    eta = sigma^2;
    h_N = exp(-abs(x).^2/eta) / eta^N;
    
end

