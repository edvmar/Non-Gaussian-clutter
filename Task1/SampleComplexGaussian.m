%%%%%%%%%%%%%% SampleComplexGaussian %%%%%%%%%%%%%%%%
%
% Returns sample from Complex Gaussian distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sample = SampleComplexGaussian(sampleSize, mean, sigma)
    sigma1dim = sigma/sqrt(2); 
    

    samples = randn(1, sampleSize); 
    a = samples*sigma1dim;
    
    samples = randn(1, sampleSize); 
    b = samples*sigma1dim;
    
    sample = a + 1i*b + mean;

end

