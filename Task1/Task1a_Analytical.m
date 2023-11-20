%%%%%%%%%%%%%%%%%%% Task 1 a Analytical %%%%%%%%%%%%%%%%%%%%
%
% Produces ROC curves for the 0D - problem analytically
% Gaussian detector and Gaussian clutter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc


SIRs = [0, 3, 10, 13]; % dB 

numberOfEtaValues = 50;

etaValues = {linspace(0.5, 10, numberOfEtaValues), linspace(0.5, 10, numberOfEtaValues),...
    linspace(0.5, 10, numberOfEtaValues), linspace(0.5, 10, numberOfEtaValues)}; % check later! 


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

    for iEta=1:numberOfEtaValues
        eta = etaValues{iSIR}(iEta); 
        a_l = (log(eta)+alpha)^2/(2*alpha);

        pFalseAlarm(iSIR,iEta) = 1 - normcdf(a_l);
        pDetection(iSIR,iEta) = 1 - normcdf(2*(a_l-alpha));
        
    end
end 

%% Plotting 
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
end
set(gca, 'XScale', 'log');
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location='northwest')

