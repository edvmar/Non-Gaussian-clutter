function h_N = TailDistributionCompoundGaussian(x, N, sigma, nu)
    
    % Tail distribution with derivatives
    eta = sigma^2;
    h_N = (nu/eta)^N/gamma(nu)*2*(sqrt(nu/eta).*abs(x)).^(nu-N)...
          .*besselk(nu-N,2*sqrt(nu/eta)*abs(x));  

end