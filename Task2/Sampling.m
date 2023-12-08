%%%%%%%%%%%%%% SampleComplexGaussionRow %%%%%%%%%%%%%%%%%
%
% Numerically samples a row of the CPI from CN(0,sigma^2)
% rMax is the largest radii that we calculate the inverse for
% Calculate the inverse numerically and choose closest value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rangeBin = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F)
    
    
    % Table for cdf x,y values
    domain = linspace(0, rMax, rMax*10);
    range = F(domain)';
    range_rep = repmat(range,1,sampleSize)';

    uniformSample = rand(numberOfPulses, sampleSize);
    xSample = zeros(numberOfPulses, sampleSize);
    

    % Sample x
    for j = 1:numberOfPulses
        uniform = uniformSample(j,:)';
        [~, index] = min(abs((range_rep-uniform)'));
        radius = domain(index);
        
        thetas = rand(1,sampleSize)*2*pi;
        xSample(j,:) = exp(thetas*1i).*radius; 
    end

    % d_k = L x
    rangeBin = (L*xSample);

end

