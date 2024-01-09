%%%%%%%%%%%%%%%%%%% Task 1 c %%%%%%%%%%%%%%%%%%%%
%
% Produces ROC curves for the 0D - problem 
% CG detector and CN clutter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc


SIRs = [0, 3, 10, 13]; % dB 

numberOfEtaValues = 1000;
maxEta = 1e4;
etaValues = {linspace(0.5, 500, numberOfEtaValues),...
             [linspace(0.5, 200, numberOfEtaValues*0.5),linspace(200, 1e4, numberOfEtaValues*0.5)],...
             [linspace(0.5, 200, numberOfEtaValues*0.3),linspace(200, maxEta, numberOfEtaValues*0.7)], ...
             [linspace(0.5, 200, numberOfEtaValues*0.3),linspace(200, maxEta, numberOfEtaValues*0.7)]}; 


sampleSize = 10^8;

detectorSigma = 1; % The standard deviation for the detector
clutterSigma = 1; % The standard deviation for the detector
detectorMean = 0;
clutterMean = 0;

sumFA = zeros(length(SIRs), numberOfEtaValues);
sumTD = zeros(length(SIRs), numberOfEtaValues);

tic
for iSIR = 1:length(SIRs)
    iSIR
    SIR = 10^(SIRs(iSIR)/10);           
    alpha = clutterSigma*sqrt(SIR);             
   
    theta = 0;
    s = alpha*(cos(theta)+1i*sin(theta)); % signal 

    clutterSample = SampleComplexGaussian(sampleSize, clutterMean, clutterSigma); 
    signalSample = clutterSample + s;

    % False Alarm (*)
    fH1_fa = CompoundGaussianPDF(clutterSample, detectorMean + s, detectorSigma);      
    fH0_fa = CompoundGaussianPDF(clutterSample, detectorMean, detectorSigma);
    LRT_fa = fH1_fa./fH0_fa;

    % True Detection (**)
    fH1_td = CompoundGaussianPDF(signalSample, detectorMean + s, detectorSigma);   
    fH0_td = CompoundGaussianPDF(signalSample, detectorMean, detectorSigma);
    LRT_td = fH1_td./fH0_td;

    for iEta=1:numberOfEtaValues
        eta = etaValues{iSIR}(iEta); 

        % False Alarm (*)
        sumFA(iSIR, iEta) = sum((LRT_fa > eta));
        
        % True Detection (**)
        sumTD(iSIR, iEta) = sum((LRT_td > eta));

    end
end 
toc

pFalseAlarm = sumFA/sampleSize;
pDetection  = sumTD/sampleSize;

%% Plotting 
figure(2)
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
end
set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location = 'west', FontSize=14)
axis([1e-7, 1, 0, 1])
 

