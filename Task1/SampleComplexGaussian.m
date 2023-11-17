%%%%%%%%%%%%%% SampleComplexGaussian %%%%%%%%%%%%%%%%
%
% Returns sample from Complex Gaussian distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sample = SampleComplexGaussian(sampleSize, mean, sigma)
    
    samples = randn(1, sampleSize); 
    a = samples*sigma + mean;

    samples = randn(1, sampleSize); 
    b = samples*sigma + mean;

    sample = a + 1i*b;

end

