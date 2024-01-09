%%%%%%%%%%%%%%%%%%%%%%%%%%% Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Samples a row/rangebin of the CPI using the numerical inverse
% of the CDF (F).
% rMax is the largest radii that we calculate the inverse for.
% L is the Cholesky decomposition of the covariance matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sampleCPI = Sampling3D(numberOfPulses, numberOfDistances, sampleSize, rMax, L, F)
    
    % Table for cdf x,y values
    domain = linspace(0, rMax, 7);
    range = F(domain)';
    rangeMatrix = repmat(range, 1, numberOfPulses, numberOfDistances, sampleSize);
    rangeMatrix = permute(rangeMatrix, [2 3 4 1]); % to get correct dimensions to calc rangeDifference
    
    % sampling
    uniformSample = rand(numberOfPulses, numberOfDistances, sampleSize);
    
    % Find inverse from table
    rangeDifference = rangeMatrix-uniformSample;
    rangeDifference = permute(rangeDifference, [4, 1, 2, 3]); % to get correct dimensions for the min operation
    [~, index] = min(abs(rangeDifference));
    radius = domain(index);
    radius = squeeze(radius); % removes the redundant first dimension of radius
    
    % rewrite sample as a complex vectors
    xSample = exp(rand(numberOfPulses, numberOfDistances, sampleSize)*2*pi*1i).*radius;
    % thetas  = rand(numberOfPulses, numberOfDistances, sampleSize)*2*pi; 
    % xSample = exp(thetas*1i).*radius; % OBS! memory expensive to do it with two steps and thus commented
    
    % Calculate the rangeBins (but since we calculate all of them at once
    % we get all the CPI samples at once)
    sampleCPI = pagemtimes(L, xSample); % pagewise matrix multiplication
                                        % each page corresponds to a CPI
                                        % matrix sample
                                        % OBS! The CPI matrices are the
                                        % transponate of the definition (in the
                                        % other codes we take the transponate
                                        % after the sampling function in the code)

end

