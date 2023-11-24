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
etaValues = {linspace(0.5, 1000, numberOfEtaValues),...
             [linspace(0.5, 200, numberOfEtaValues*0.1),linspace(200, 1e4, numberOfEtaValues*0.9)],...
             [linspace(0.5, 1000, numberOfEtaValues*0.5),linspace(1000, maxEta, numberOfEtaValues*0.5)], ...
             [linspace(0.5, 1000, numberOfEtaValues*0.5),linspace(1000, maxEta, numberOfEtaValues*0.5)]}; 
    


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

%% P_FA vs eta 
figure(2)
subplot(1,3,1)
hold on

plot(etaValues{1},pFalseAlarm(1,:), LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('SIR = 0', location='best')
axis([0, 1e3, 1e-8, 1])

subplot(1,3,2)
hold on
set(gca,'ColorOrderIndex',2)
plot(etaValues{2},pFalseAlarm(2,:), LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('SIR = 3', location='best')
axis([0, 1e4, 1e-8, 1])

subplot(1,3,3)
hold on
set(gca,'ColorOrderIndex',3)
for iSIR = 3:4
    plot(etaValues{iSIR},pFalseAlarm(iSIR,:), LineWidth=1.5)
end
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('SIR = 10', 'SIR = 13', location='best')
axis([0, maxEta, 1e-8, 1])









