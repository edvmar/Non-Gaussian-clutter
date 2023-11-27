%%%%%%%%%%%%%%%%%%% Task 1 a %%%%%%%%%%%%%%%%%%%%
%
% Produces ROC curves for the 0D - problem 
% Gaussian detector and Gaussian clutter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc


SIRs = [0, 3, 10, 13]; % dB 

numberOfEtaValues = 1000;
maxEta = 5*1e5;
etaValues = {linspace(0.5, 1000, numberOfEtaValues),...
             [linspace(0.5, 200, numberOfEtaValues*0.1),linspace(200, 1e4, numberOfEtaValues*0.9)],...
             [linspace(0.5, 1000, numberOfEtaValues*0.5),linspace(1000, maxEta, numberOfEtaValues*0.5)], ...
             [linspace(0.5, 1000, numberOfEtaValues*0.5),linspace(1000, maxEta, numberOfEtaValues*0.5)]}; 


sampleSize = 10^4; % 10^8 later? 

detectorSigma = 1; % The standard deviation for the detector
clutterSigma = 1; % The standard deviation for the detector
detectorMean = 0;
clutterMean = 0;

% pFalseAlarm = zeros(length(SIRs), numberOfEtaValues);
% pDetection = zeros(length(SIRs), numberOfEtaValues);
sumFA = zeros(length(SIRs), numberOfEtaValues);
sumTD = zeros(length(SIRs), numberOfEtaValues);

tic
for iSIR = 1:length(SIRs)
    
    SIR = 10^(SIRs(iSIR)/10);           
    alpha = clutterSigma*sqrt(SIR);             
   
    theta = 0; % change to rand(1,1)*2*pi ? 
    s = alpha*(cos(theta)+1i*sin(theta)); % signal 

    clutterSample = SampleComplexGaussian(sampleSize, clutterMean, clutterSigma); 
    signalSample = clutterSample + s;

    % False Alarm (*)
    fH1_fa = ComplexGaussianPDF(clutterSample, detectorMean + s, detectorSigma);           % or clutter mean?
    fH0_fa = ComplexGaussianPDF(clutterSample, detectorMean, detectorSigma);

    % True Detection (**)
    fH1_td = ComplexGaussianPDF(signalSample, detectorMean + s, detectorSigma);           % or clutter mean?
    fH0_td = ComplexGaussianPDF(signalSample, detectorMean, detectorSigma);

    for iEta=1:numberOfEtaValues
        eta = etaValues{iSIR}(iEta); 

        % False Alarm (*)
        sumFA(iSIR, iEta) = sum(((fH1_fa./fH0_fa) > eta));
        % sumFA = sum(((fH1_fa./fH0_fa) > eta));
        
        % True Detection (**)
        sumTD(iSIR, iEta) = sum(((fH1_td./fH0_td) > eta));
        %sumTD = sum(((fH1_td./fH0_td) > eta));


        % pFalseAlarm(iSIR, iEta) = sumFA/sampleSize;
        % pDetection(iSIR, iEta) = sumTD/sampleSize;
        iEta
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
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location = 'southeast')
axis([1e-7, 1, 0, 1])
 

