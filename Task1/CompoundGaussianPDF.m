%%%%%%%%%%%%%% CompoundGaussianPDG %%%%%%%%%%%%%%%%
%
% Returns pdf values for Complex Gauss distribution
% given mean value and shape parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = CompoundGaussianPDF(z, mean, sigma)
    
    %b = shape / mean;

    %f = 2*b./(gamma(shape)*pi).* (sqrt(b)*abs(z)).^(shape-1).*besselk(shape-1, 2*sqrt(b)*abs(z));

    % REWRITE WITH ANALYTICAL FINDING:

    % f = besselk(0, 2*abs(z-mean)/sigma)/(2*pi*(sigma/2)^2);


    % Weibull
    f = exp(-abs(z-mean)*sqrt(2)/sigma)./(pi*sqrt(2)*sigma*abs(z-mean));

end

