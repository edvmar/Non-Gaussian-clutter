function h_N = H_nKdist(y, numberOfPulses, sigma, nu)
    
    % Tail distribution with derivatives
    N = numberOfPulses;
    eta = sigma^2;

    h_N = (nu/eta)^N/gamma(nu)*2*sqrt(nu/eta*y).^(nu-N)...
          .*besselk(nu-N,2*sqrt(nu/eta*y));  
    
    % Blir inte det här wack när N > nu ? N - nu i modern, N - nu i
    % Sangston 

    %b = nu/eta;
    %h_N = (2./y).*(b*y).^((N+nu)/2)./(gamma(N)*gamma(nu)).*besselk(nu-N, 2*sqrt(b*y));

end