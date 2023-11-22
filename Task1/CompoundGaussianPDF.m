%%%%%%%%%%%%%% CompoundGaussianPDG %%%%%%%%%%%%%%%%
%
% Returns pdf values for Complex Gauss distribution
% given mean value and shape parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = CompoundGaussianPDF(z, mean, shape)
    
    b = shape / mean;

    f = 2*b./(gamma(shape)*pi).* (sqrt(b)*abs(z)).^(shape-1).*besselk(shape-1, 2*sqrt(b)*abs(z));

    % REWRITE WITH ANALYTICAL FINDING

end

