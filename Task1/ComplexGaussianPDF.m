%%%%%%%%%%%%%% ComplexGaussianPDF %%%%%%%%%%%%%%%%
%
% Returns pdf values for Complex Gaussian distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = ComplexGaussianPDF(z, mean, sigma)

    f = 1/(pi*sigma^2).*exp(-abs(z-mean).^2/sigma^2); 

end

