%%%%%%%%%%%%%% SampleCompoundGaussian %%%%%%%%%%%%%%%%%%%
%
% Returns sample from the Compound Gaussian distribution given by
% Z \sim CN(0,S) where S \sim N(0,\sigma)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sample = SampleCompoundGaussian(sampleSize, mean, sigma)
    
    sigmaSample = abs(randn(1,sampleSize)*sigma);
    
    a = randn(1,sampleSize).*sigmaSample/sqrt(2);
    b = randn(1,sampleSize).*sigmaSample/sqrt(2);
   
    sample = a + 1i*b + mean; 

end