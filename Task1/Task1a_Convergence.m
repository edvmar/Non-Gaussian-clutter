%%%%%%%%%%%%%%%%%%% Task 1 a Convergence %%%%%%%%%%%%%%%%%%%%
%
% Analyzing the convergence of P_FA with increasing sampleSize
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

sampleSizes = fix(logspace(3,8,200));

eta = 100;

detectorSigma = 1; % The standard deviation for the detector
clutterSigma  = 1; % The standard deviation for the detector
detectorMean  = 0;
clutterMean   = 0;

SIR = 10; % dB 
SIR = 10^(SIR/10); 
alpha = clutterSigma*sqrt(SIR);             

theta = 0; 
s = alpha*(cos(theta)+1i*sin(theta)); % signal 

pFalseAlarmAnalytical = zeros(length(sampleSizes), 1);
pFalseAlarmNumerical = zeros(length(sampleSizes), 1);


for iSampleSize = 1:length(sampleSizes)
    iSampleSize

    sampleSize = sampleSizes(iSampleSize);

    clutterSample = SampleComplexGaussian(sampleSize, clutterMean, clutterSigma); 
    signalSample = clutterSample + s;
    
    % Numerical
    fH1_fa = ComplexGaussianPDF(clutterSample, detectorMean + s, detectorSigma);        
    fH0_fa = ComplexGaussianPDF(clutterSample, detectorMean, detectorSigma);
    sumFA = sum(((fH1_fa./fH0_fa) > eta));

    pFalseAlarmNumerical(iSampleSize) = sumFA/sampleSize;

    % Analytical 
    threshold = (log(eta)+alpha^2)/(2*alpha);
    pFalseAlarmAnalytical(iSampleSize) = 1 - normcdf(sqrt(2)*threshold);

end 

errorFalseAlarm = abs(pFalseAlarmAnalytical - pFalseAlarmNumerical);

%% Plotting 
figure(1)
hold on
plot(sampleSizes, errorFalseAlarm, LineWidth=1.5)
plot(sampleSizes,1./sqrt(sampleSizes), LineWidth=1.5)
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Error $P_{FA}$', '$\frac{1}{\sqrt{M}}$', 'Interpreter', 'latex', FontSize=14) 
xlabel('M'), ylabel('Error')

