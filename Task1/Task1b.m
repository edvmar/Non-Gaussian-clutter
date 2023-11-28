%%%%%%%%%%%%%%%%%%% Task 1 b %%%%%%%%%%%%%%%%%%%%
%
% Produces ROC curves for the 0D - problem 
% Gaussian detector and Compound clutter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc


SIRs = [0, 3, 10, 13]; % dB 
numberOfEtaValues = 2000;
secLast = 0;
probSecLast = [0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.03, 0.05, 0.06, 0.09, 0.11, 0.14, 0.15, 0.15, 0.14] ;
last = 0;
probLast = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.04, 0.04,...
    0.05, 0.05, 0.05, 0.06, 0.06, 0.08, 0.09, 0.09, 0.09];
for i=1:length(probSecLast)
    secLast = [secLast, linspace(10^(i-1), 10^i, numberOfEtaValues*probSecLast(i))];
end

for i=1:length(probLast)
    last = [last, linspace(10^(i-1), 10^i, numberOfEtaValues*probLast(i))];
end

secLast = secLast + 0.5;
last = last + 0.5;
etaValues = {[linspace(0.5, 100, numberOfEtaValues*0.8),linspace(100, 1e4, numberOfEtaValues*0.2)],...
             [linspace(0.5, 100, numberOfEtaValues*0.1),linspace(100, 1e4, numberOfEtaValues*0.9)],...
             secLast, ...
             last}; 


sampleSize = 10^4; % 10^8 later? 

detectorSigma = 1; % The standard deviation for the detector
clutterSigma = 1; % The standard deviation for the detector
detectorMean = 0;
clutterMean = 0;

% pFalseAlarm = zeros(length(SIRs), numberOfEtaValues);
% pDetection = zeros(length(SIRs), numberOfEtaValues);
sumFA = zeros(length(SIRs), numberOfEtaValues);
sumTD = zeros(length(SIRs), numberOfEtaValues);
iter = 0;
tic
for iSIR = 1:length(SIRs)
    iSIR
    SIR = 10^(SIRs(iSIR)/10);           
    alpha = clutterSigma*sqrt(SIR);             
   
    theta = 0; % change to rand(1,1)*2*pi ? 
    s = alpha*(cos(theta)+1i*sin(theta)); % signal 

    clutterSample = SampleCompoundGaussian(sampleSize, clutterMean, clutterSigma); 
    signalSample = clutterSample + s;

    % False Alarm (*)
    fH1_fa = ComplexGaussianPDF(clutterSample, detectorMean + s, detectorSigma);           % or clutter mean?
    fH0_fa = ComplexGaussianPDF(clutterSample, detectorMean, detectorSigma);
    LRT_fa = fH1_fa./fH0_fa;

    % True Detection (**)
    fH1_td = ComplexGaussianPDF(signalSample, detectorMean + s, detectorSigma);           % or clutter mean?
    fH0_td = ComplexGaussianPDF(signalSample, detectorMean, detectorSigma);
    LRT_td = fH1_td./fH0_td;


    for iEta=1:numberOfEtaValues
        eta = etaValues{iSIR}(iEta); 
        
        % False Alarm (*)
        sumFA(iSIR, iEta) = sum((LRT_fa > eta));
        
        % True Detection (**)
        sumTD(iSIR, iEta) = sum((LRT_td > eta));

    end
    
end 
pFalseAlarm = sumFA/sampleSize;
pDetection  = sumTD/sampleSize;
toc

%% Plotting 
figure(2)
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
end
set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location = 'west')
axis([1e-7, 1, 0, 1])
 

