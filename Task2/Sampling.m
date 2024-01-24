%%%%%%%%%%%%%%%%%%%%%%%%%%% Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Samples a row/rangebin of the CPI using the numerical inverse
% of the CDF (F).
% rMax is the largest radii that we calculate the inverse for.
% L is the Cholesky decomposition of the covariance matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rangeBin = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F)
    
    % Table for cdf x,y values
    domain = linspace(0, rMax, rMax*100);
    range = F(domain)';
    range_rep = repmat(range,1,sampleSize)'; % NOTE: Very data-inefficient. See ReadMe for suggested change.

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

    % d_k = L x to get covariance
    rangeBin = (L*xSample);

end

