%%%%%%%%%%%%%%%%%%% Task 1 a %%%%%%%%%%%%%%%%%%%%
%
% Produces ROC curves for the 0D - problem 
% Gaussian detector and Gaussian clutter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc


SIRs = [0, 3, 10, 13]; % dB 

numberOfEtaValues = 10;

etaValues = {linspace(1, 10, numberOfEtaValues), linspace(1, 10, numberOfEtaValues),...
    linspace(1, 10, numberOfEtaValues), linspace(1, 10, numberOfEtaValues)}; % check later! 


sampleSize = 10^6; % 10^8 later? 

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
    alpha = clutterSigma*sqrt(SIR);             
   
    theta = 0;
    s = alpha*(cos(theta)+1i*sin(theta)); % signal 

    clutterSample = SampleComplexGaussian(sampleSize, clutterMean, clutterSigma); % ?
    signalSample = clutterSample + s;

    for iEta=1:numberOfEtaValues
        eta = etaValues{iSIR}(iEta); 

        
        % False Alarm 
        fH1_fa = ComplexGaussianPDF(clutterSample, detectorMean + s, detectorSigma);           % or clutter mean?
        fH0_fa = ComplexGaussianPDF(clutterSample, detectorMean, detectorSigma);
        sumFA = sum(((fH1_fa./fH0_fa) > eta));
        
        % True Detection
        fH1_td = ComplexGaussianPDF(signalSample, detectorMean + s, detectorSigma);           % or clutter mean?
        fH0_td = ComplexGaussianPDF(signalSample, detectorMean, detectorSigma);
        sumTD = sum(((fH1_td./fH0_td) > eta));

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
%% Analytically
clc,clear
hold on
numberOfEtaValues = 20;
etaValues = linspace(1, 10, numberOfEtaValues);
SIRs = 8;
clutterSigma=1;
pFA = zeros(1, numberOfEtaValues);
pTD = zeros(1, numberOfEtaValues);
for iEta=1:numberOfEtaValues
    eta = etaValues(iEta); 
    SIR = 10^(SIRs/10);   
    alpha = clutterSigma*sqrt(SIR);             
    a_l = (log(eta)+alpha)^2/(2*alpha);
    pFA(iEta) = 1 - normcdf(a_l);
    pTD(iEta) = 1 - normcdf(2*(a_l-alpha));
end

plot(pFA, pTD, LineWidth=1.5)
set(gca, 'XScale', 'log');
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13')


