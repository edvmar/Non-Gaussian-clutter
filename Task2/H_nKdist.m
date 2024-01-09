%%%%%%%%%%%%%%%%%%%% H_nGaussian %%%%%%%%%%%%%%%%%%%%%%%
%
% Tail distribution for the complex K distribution 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h_N = H_nKdist(y, numberOfPulses, sigma, nu)
    
    N = numberOfPulses;
    eta = sigma^2;

    h_N = (nu/eta)^N/gamma(nu)*2*sqrt(nu/eta*y).^(nu-N)...
          .*besselk(nu-N,2*sqrt(nu/eta*y));  

end