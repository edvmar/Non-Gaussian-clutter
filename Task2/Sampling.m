%%%%%%%%%%%%%% SampleComplexGaussionRow %%%%%%%%%%%%%%%%%
%
% Numerically samples a row of the CPI from CN(0,sigma^2)
% rMax is the largest radii that we calculate the inverse for
% Calculate the inverse numerically and choose closest value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rangeBin = Sampling(numberOfPulses, rMax, L, F)
    
    
    % Table for cdf x,y values
    domain = linspace(0, rMax, rMax*1000);
    range = F(domain);

    uniformSample = rand(1, numberOfPulses);
    xSample = zeros(numberOfPulses, 1);
    
    % Sample x
    for j = 1:numberOfPulses
        uniform = uniformSample(j);
        [~, index] = min(abs(range-uniform));
        radius = domain(index);
        
        thetas = rand(1)*2*pi;
        xSample(j) = exp(thetas*1i)*radius; 
    end

    % d_k = L x
    rangeBin = (L*xSample)';

end

