%%%%%%%%%%%%%%%%%%% Task 1 a %%%%%%%%%%%%%%%%%%%%
%
% Produces ROC curves for the 0D - problem 
% Gaussian detector and Gaussian clutter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

SIRs = [1, 3, 5, 8]; % dB 

numberOfEtaValues = 40;
etaValues = {linspace(0.01, 10, numberOfEtaValues),linspace(0.01, 10, numberOfEtaValues), 
    linspace(0.01, 10, numberOfEtaValues), linspace(0.01, 10, numberOfEtaValues)}; % check later! 

sampleSize = 10^5; % 10^8 later? 

detectorSigma = 1; % The standard deviation for the detector
clutterSigma = 1; % The standard deviation for the detector
detectorMean = 0;
clutterMean = 0;

pFalseAlarm = zeros(length(SIRs), numberOfEtaValues);
pDetection = zeros(length(SIRs), numberOfEtaValues);

for iSIR = 1:length(SIRs)
    % SIR = SIRs(iSIR); 
    % alpha = clutterSigma*sqrt(SIR);    % signal strength  dunno if this is correct?
    SIR = 10^(SIRs(iSIR)/10);           % potentially like this ? 
    alpha = clutterSigma*sqrt(SIR);              % Samuel question
   
    theta = 0;
    s = alpha*(cos(theta)+1i*sin(theta)) % signal 

    clutterSample = SampleComplexGaussian(sampleSize, clutterMean, clutterSigma);
    signalSample = clutterSample + s;

    for iEta=1:numberOfEtaValues
        eta = etaValues{iSIR}(iEta); 
        sumFA = 0;
        sumTD = 0;

        for j=1:sampleSize
            
            % False Alarm
            fH1 = ComplexGaussianPDF(clutterSample(j), detectorMean + s, detectorSigma);           % or clutter mean?
            fH0 = ComplexGaussianPDF(clutterSample(j), detectorMean, detectorSigma);
            if (fH1/fH0 > eta)
                sumFA = sumFA + 1;
            end
            
            % True Detection
            fH1 = ComplexGaussianPDF(signalSample(j), detectorMean + s, detectorSigma);           % or clutter mean?
            fH0 = ComplexGaussianPDF(signalSample(j), detectorMean, detectorSigma);
            if (fH1/fH0 > eta)
                sumTD = sumTD + 1;
            end
        end

        pFalseAlarm(iSIR, iEta) = sumFA/sampleSize;
        pDetection(iSIR, iEta) = sumTD/sampleSize;
    end

end 

%% Plotting 
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
end
set(gca, 'XScale', 'log');
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13')

%Should probably be with log axis etc.. 
