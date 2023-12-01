function clutterRow = SampleComplexGaussianRow(numberOfPulses, rMax, sigma, L)
    
    F = @(x) 1 - exp(-abs(x).^2);  % eqn (12) sigma? 

    xValues = linspace(0, rMax, rMax*1000);
    yValues = F(xValues);

    uniformSample = rand(1,numberOfPulses);
    xSample = zeros(numberOfPulses,1);
    
    % Sample x
    
    for j = 1:numberOfPulses
        uniform = uniformSample(j);
        [~, index] = min(abs(yValues-uniform));
        radius = xValues(index);
        
        thetas = rand(1)*2*pi;
        xSample(j) = exp(thetas*1i)*radius; 
    end

    clutterRow = L*xSample;

end

