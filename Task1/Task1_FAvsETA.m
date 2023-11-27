%%%%%%%%%%%%%%%%%%% Task 1 P_FA and Eta %%%%%%%%%%%%%%%%%%%%
%
% Analyzing the relation between P_FA and eta for the four cases
% with a specified SIR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

sampleSize = 2*1e6;

numberOfEtaValues = 1000;
maxEta = 4*1e5;
etaValues = [linspace(0.5, 5000, numberOfEtaValues*0.5),linspace(5000, maxEta, numberOfEtaValues*0.5)]; 


detectorSigma = 1; % The standard deviation for the detector
clutterSigma = 1; % The standard deviation for the detector
detectorMean = 0;
clutterMean = 0;


SIR = 10; % dB
SIR = 10^(SIR/10);         
alpha = clutterSigma*sqrt(SIR);
theta = 0; 
s = alpha*(cos(theta)+1i*sin(theta)); % signal 


clutterSampleCN = SampleComplexGaussian(sampleSize, clutterMean, clutterSigma); % ?
%signalSampleCN = clutterSampleCN + s;

clutterSampleCG = SampleCompoundGaussian(sampleSize, clutterMean, clutterSigma); % ?
%signalSampleCG = clutterSampleCG + s;



pFalseAlarmCNCN = zeros(1, numberOfEtaValues);
pFalseAlarmCNCG = zeros(1, numberOfEtaValues);
pFalseAlarmCGCN = zeros(1, numberOfEtaValues);
pFalseAlarmCGCG = zeros(1, numberOfEtaValues);


for iEta=1:numberOfEtaValues
    eta = etaValues(iEta); 
    
    % False Alarm 
    fH1_fa = ComplexGaussianPDF(clutterSampleCN, detectorMean + s, detectorSigma); 
    fH0_fa = ComplexGaussianPDF(clutterSampleCN, detectorMean, detectorSigma);
    sumFA_CNCN = sum(((fH1_fa./fH0_fa) > eta));

    % False Alarm 
    fH1_fa = ComplexGaussianPDF(clutterSampleCG, detectorMean + s, detectorSigma); 
    fH0_fa = ComplexGaussianPDF(clutterSampleCG, detectorMean, detectorSigma);
    sumFA_CNCG = sum(((fH1_fa./fH0_fa) > eta));

    % False Alarm 
    fH1_fa = CompoundGaussianPDF(clutterSampleCN, detectorMean + s, detectorSigma); 
    fH0_fa = CompoundGaussianPDF(clutterSampleCN, detectorMean, detectorSigma);
    sumFA_CGCN = sum(((fH1_fa./fH0_fa) > eta));

    % False Alarm 
    fH1_fa = CompoundGaussianPDF(clutterSampleCG, detectorMean + s, detectorSigma); 
    fH0_fa = CompoundGaussianPDF(clutterSampleCG, detectorMean, detectorSigma);
    sumFA_CGCG = sum(((fH1_fa./fH0_fa) > eta));
    


    pFalseAlarmCNCN(iEta) = sumFA_CNCN/sampleSize;
    pFalseAlarmCNCG(iEta) = sumFA_CNCG/sampleSize;
    pFalseAlarmCGCN(iEta) = sumFA_CGCN/sampleSize;
    pFalseAlarmCGCG(iEta) = sumFA_CGCG/sampleSize;

end


%% Plotting
subplot(1,2,1)
hold on
plot(etaValues,pFalseAlarmCNCN, LineWidth=1.5)
plot(etaValues,pFalseAlarmCNCG, LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('CN-CN', 'CN-CG', location='best')
axis([0, maxEta, 1e-8, 1])

subplot(1,2,2)
set(gca,'ColorOrderIndex',3)
hold on
plot(etaValues,pFalseAlarmCGCN, LineWidth=1.5)
plot(etaValues,pFalseAlarmCGCG, LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('CG-CN','CG-CG', location='best')
axis([0, 1e4, 1e-8, 1])




