function h_N = TailDistributionCompoundGaussian(y, numberOfPulses, sigma, nu)
    
    % Tail distribution with derivatives
    N = numberOfPulses;
    eta = sigma^2;
    h_N = (nu/eta)^N/gamma(nu)*2*sqrt(nu/eta*y).^(nu-N)...
          .*besselk(nu-N,2*sqrt(nu/eta*y));  
    
    % Blir inte det här wack när N > nu ? N - nu i modern, N - nu i
    % Sangston
end