%%%%%%%%%%%%%% CompoundGaussianPDG %%%%%%%%%%%%%%%%
%
% Returns pdf values for Complex Gauss distribution
% given mean value and standard deviation.
%
% Here a Weibull distribution was derived as the distribution of 
% Z \sim CN(0,S) where S \sim N(0,\sigma)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = CompoundGaussianPDF(z, mean, sigma)
    
    f = exp(-abs(z-mean)*sqrt(2)/sigma)./(pi*sqrt(2)*sigma*abs(z-mean));

end

