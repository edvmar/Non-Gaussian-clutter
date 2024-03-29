%%%%%%%%%%%%%%%%%%% Task 1 a Numerical and Analytical %%%%%%%%%%
%
% Numerical and analytical ROC curves for the CN-CN case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

sampleSize = 1e8;

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

pFalseAlarmAnalytical = zeros(length(SIRs), numberOfEtaValues);
pDetectionAnalytical = zeros(length(SIRs), numberOfEtaValues);

pFalseAlarmNumerical = zeros(length(SIRs), numberOfEtaValues);
pDetectionNumerical = zeros(length(SIRs), numberOfEtaValues);

for iSIR = 1:length(SIRs)
    
    SIR = 10^(SIRs(iSIR)/10);      
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

        pFalseAlarmNumerical(iSIR, iEta) = sumFA/sampleSize;
        pDetectionNumerical(iSIR, iEta) = sumTD/sampleSize;
    
        % Analytical 
        threshold = (log(eta)+alpha^2)/(2*alpha);

        pFalseAlarmAnalytical(iSIR,iEta) = 1 - normcdf(sqrt(2)*threshold);
        pDetectionAnalytical(iSIR,iEta) = 1 - normcdf(sqrt(2)*(threshold-alpha));
        
    end


end 

%% Plotting 
figure(1)
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarmAnalytical(iSIR,:), pDetectionAnalytical(iSIR, :), LineWidth=3)
end
set(gca,'ColorOrderIndex',1)
for iSIR = 1:length(SIRs)
    plot(pFalseAlarmNumerical(iSIR,:), pDetectionNumerical(iSIR, :),'k--', LineWidth=2)
end
set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', 'Numerical', location='southeast')
axis([1e-7, 1, 0, 1])


%% P_FA vs eta 
figure(2)
subplot(1,3,1)
hold on

plot(etaValues{1},pFalseAlarmAnalytical(1,:), LineWidth=1.5)
set(gca,'ColorOrderIndex',1)
plot(etaValues{1},pFalseAlarmNumerical(1,:), '--', LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('SIR = 0', location='best')
axis([0, 1e3, 1e-8, 1])

subplot(1,3,2)
hold on
set(gca,'ColorOrderIndex',2)
plot(etaValues{2},pFalseAlarmAnalytical(2,:), LineWidth=1.5)
set(gca,'ColorOrderIndex',2)
plot(etaValues{2},pFalseAlarmNumerical(2,:), '--', LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('SIR = 3', location='best')
axis([0, 1e4, 1e-8, 1])

subplot(1,3,3)
hold on
set(gca,'ColorOrderIndex',3)
for iSIR = 3:4
    plot(etaValues{iSIR},pFalseAlarmAnalytical(iSIR,:), LineWidth=1.5)
end
set(gca,'ColorOrderIndex',3)
for iSIR = 3:4
    plot(etaValues{iSIR},pFalseAlarmNumerical(iSIR,:), '--', LineWidth=1.5)
end
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('SIR = 10', 'SIR = 13', location='best')
axis([0, maxEta, 1e-8, 1])





