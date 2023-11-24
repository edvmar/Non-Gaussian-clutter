%%%%%%%%%%%%%% SampleCompoundGaussian %%%%%%%%%%%%%%%%
%
% Returns sample from K - distribution
% Given mean and shape parameters
% Note: Don't know sigma yet.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sample = SampleCompoundGaussian(sampleSize, mean, sigma)
    
    
    sigmaSample = abs(randn(1,sampleSize)*sigma);
    
    a = randn(1,sampleSize).*sigmaSample/sqrt(2);
    b = randn(1,sampleSize).*sigmaSample/sqrt(2);
   
    sample = a + 1i*b + mean; 
    

end

%     GaussianData = randn(1, sampleSize);
% 
%     Q = 0.5*(erfc(GaussianData/(2^0.5)));
% 
%     GammaData = gaminv(Q,shape,mean);
% 
%     modulation = GammaData.* mean;
% 
%     r = rand(sampleSize,1);
%     kData = -modulation.* log(r');
% 
%     sampleRadii = kData.* mean;
% 
%     sampleTheta = rand(1, sampleSize)*2*pi;
% 
%     sample = sampleRadii.*exp(1i*sampleTheta); 
%     % Don't know if this is how to do it complex but the Kdist only gives positive values.. 
% 
%     % if(rand(1)<0.5)
%     % *-1
% 
%     % Guessing there's a normalisation factor missing ? 