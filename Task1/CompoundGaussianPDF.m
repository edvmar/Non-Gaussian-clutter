%%%%%%%%%%%%%% CompoundGaussianPDG %%%%%%%%%%%%%%%%
%
% Returns pdf values for Complex Gauss distribution
% given mean value and shape parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = CompoundGaussianPDF(z, mean, shape)
    
    b = shape / mean;

    f = 2*b./gamma(shape).* (sqrt(b*z)).^(shape-1).*besselk(shape-1, z);

end

