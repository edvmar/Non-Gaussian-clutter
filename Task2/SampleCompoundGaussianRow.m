%%%%%%%%%%%%%% SampleComplexGaussionRow %%%%%%%%%%%%%%%%%
%
% Numerically samples a row of the CPI from CG(0,sigma^2)
% rMax is the largest radii that we calculate the inverse for
% Calculate the inverse numerically and choose closest value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clutterRow = SampleCompoundGaussianRow(numberOfPulses, rMax, nu, sigma, L)
    
    % cdf
    eta = sigma^2;
    F = @(x) 1 - (2*(sqrt(nu/eta).*abs(x)).^nu)/gamma(nu).*besselk(nu,2*sqrt(nu/eta)*abs(x));  % eqn (12)   (maybe nu-1 or nu in Bessel ??)

    % table for cdf x,y values 
    xValues = linspace(0, rMax, rMax*1000);
    yValues = F(xValues);

    uniformSample = rand(1,numberOfPulses);
    xSample = zeros(numberOfPulses,1);
    
    % Sample x
    for j = 1:numberOfPulses
        uniform = uniformSample(j);
        [~, index] = min(abs(yValues-uniform));
        radius = xValues(index);
        
        thetas = rand(1)*2*pi;
        xSample(j) = exp(thetas*1i)*radius; 
    end

    % d_k = L x
    clutterRow = L*xSample;

end

