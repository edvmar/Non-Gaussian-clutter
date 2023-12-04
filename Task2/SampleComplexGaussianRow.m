%%%%%%%%%%%%%% SampleComplexGaussionRow %%%%%%%%%%%%%%%%%
%
% Numerically samples a row of the CPI from CN(0,sigma^2)
% rMax is the largest radii that we calculate the inverse for
% Calculate the inverse numerically and choose closest value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rangeBin = SampleComplexGaussianRow(numberOfPulses, rMax, sigma, L)
    
    % cdf
    F = @(x) 1 - TailDistributionComplexGaussian(x, 0, sigma);  % eqn (12) 
    
    % Table for cdf x,y values
    xValues = linspace(0, rMax, rMax*1000);
    yValues = F(xValues);

    uniformSample = rand(1, numberOfPulses);
    xSample = zeros(numberOfPulses, 1);
    
    % Sample x
    for j = 1:numberOfPulses
        uniform = uniformSample(j);
        [~, index] = min(abs(yValues-uniform));
        radius = xValues(index);
        
        thetas = rand(1)*2*pi;
        xSample(j) = exp(thetas*1i)*radius; 
    end

    % d_k = L x
    rangeBin = L*xSample;

end

