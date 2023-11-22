%%%%%%%%%%%%%%%%%%% Task 1 a Analytical %%%%%%%%%%%%%%%%%%%%
%
% Produces ROC curves for the 0D - problem analytically
% Gaussian detector and Gaussian clutter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc


SIRs = [0, 3, 10, 13]; % dB 

numberOfEtaValues = 1000;

maxEta = 4*1e5;
etaValues = {[linspace(0.5, 200, numberOfEtaValues*0.5),linspace(200, maxEta, numberOfEtaValues*0.5)]...
             [linspace(0.5, 200, numberOfEtaValues*0.5),linspace(200, maxEta, numberOfEtaValues*0.5)],...
             [linspace(0.5, 200, numberOfEtaValues*0.5),linspace(200, maxEta, numberOfEtaValues*0.5)], ...
             [linspace(0.5, 200, numberOfEtaValues*0.5),linspace(200, maxEta, numberOfEtaValues*0.5)]}; 
    


detectorSigma = 1; % The standard deviation for the detector
clutterSigma = 1; % The standard deviation for the detector
detectorMean = 0;
clutterMean = 0;

pFalseAlarm = zeros(length(SIRs), numberOfEtaValues);
pDetection = zeros(length(SIRs), numberOfEtaValues);

for iSIR = 1:length(SIRs)
    SIR = 10^(SIRs(iSIR)/10);          
    alpha = clutterSigma*sqrt(SIR);             

    for iEta=1:numberOfEtaValues
        eta = etaValues{iSIR}(iEta); 
        threshold = (log(eta)+alpha^2)/(2*alpha);

        pFalseAlarm(iSIR,iEta) = 1 - normcdf(sqrt(2)*threshold);
        pDetection(iSIR,iEta) = 1 - normcdf(sqrt(2)*(threshold-alpha));
        
    end
end 

%% Plotting 
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
end
set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location='best')
axis([1e-7, 1, 0, 1])
