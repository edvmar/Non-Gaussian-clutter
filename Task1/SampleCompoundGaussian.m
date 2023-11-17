%%%%%%%%%%%%%% SampleCompoundGaussian %%%%%%%%%%%%%%%%
%
% Returns sample from K - distribution
% Given mean and shape parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sample = SampleCompoundGaussian(sampleSize, mean, shape)
    
    samples = zeros(2, sampleSize);

    for i = 1:2
        GaussianData = randn(1, sampleSize);
    
        Q = 0.5*(erfc(GaussianData/(2^0.5)));
    
        GammaData = gaminv(Q,shape,mean);
    
        modulation = GammaData.* mean;
    
        r = rand(sampleSize,1);
        kData = -modulation.* log(r');
    
        samples(i,:) = kData.* mean;
    end

    sample = samples(1,:) + 1i*samples(2,:);     % Don't know if this is how to do it complex but.. 
end

