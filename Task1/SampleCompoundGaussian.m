%%%%%%%%%%%%%% SampleCompoundGaussian %%%%%%%%%%%%%%%%
%
% Returns sample from K - distribution
% Given mean and shape parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sample = SampleCompoundGaussian(sampleSize, mean, shape)
    
    
    GaussianData = randn(1, sampleSize);

    Q = 0.5*(erfc(GaussianData/(2^0.5)));

    GammaData = gaminv(Q,shape,mean);

    modulation = GammaData.* mean;

    r = rand(sampleSize,1);
    kData = -modulation.* log(r');

    sampleRadii = kData.* mean;

    sampleTheta = rand(1, sampleSize)*2*pi;

    sample = sampleRadii.*exp(1i*sampleTheta); 
    % Don't know if this is how to do it complex but the Kdist only gives positive values.. 


    % Guessing there's a normalisation factor missing ? 
end

