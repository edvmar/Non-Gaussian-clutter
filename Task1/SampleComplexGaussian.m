%%%%%%%%%%%%%% SampleComplexGaussian %%%%%%%%%%%%%%%%
%
% Returns sample from Complex Gaussian distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sample = SampleComplexGaussian(sampleSize, mean, sigma)
    sigma_one_dim = sigma/2; % Should it be 1/sqrt(2)?
    samples = randn(1, sampleSize); 
    a = samples*sigma_one_dim + mean;
    samples = randn(1, sampleSize); 
    b = samples*sigma_one_dim + mean;
    sample = a + 1i*b;

end

