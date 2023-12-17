%%%%%%%%%%%%%% SampleComplexGaussionRow %%%%%%%%%%%%%%%%%
%
% Numerically samples a CPI matrices from CN(0,sigma^2)
% rMax is the largest radii that we calculate the inverse for
% Calculate the inverse numerically and choose closest value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sampleCPI = Sampling3D(numberOfPulses, numberOfDistances, sampleSize, rMax, L, F)
    
    
    % create table for cdf x and y values
    domain = linspace(0, rMax, 7); % Byt till inget rMax beroende?
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
    
    % rewrite sample as a comples number
    thetas = rand(1,numberOfDistances)*2*pi;
    xSample = exp(thetas*1i).*radius; 
    
    % create CPI matrices
    sampleCPI = pagemtimes(L,xSample); % pagewise matrix multiplication each 
                                       % page corresponds to a CPI sample


end

