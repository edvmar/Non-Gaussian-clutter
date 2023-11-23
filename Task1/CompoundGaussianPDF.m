%%%%%%%%%%%%%% CompoundGaussianPDG %%%%%%%%%%%%%%%%
%
% Returns pdf values for Complex Gauss distribution
% given mean value and shape parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = CompoundGaussianPDF(z, sigma)
    
    %b = shape / mean;

    %f = 2*b./(gamma(shape)*pi).* (sqrt(b)*abs(z)).^(shape-1).*besselk(shape-1, 2*sqrt(b)*abs(z));

    % REWRITE WITH ANALYTICAL FINDING:

    f = besselk(0, abs(z)/sigma)/(2*pi*sigma^2);
    

end

