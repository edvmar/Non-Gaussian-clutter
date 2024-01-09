%%%%%%%%%%%%%%%%%%% Task 1 P_FA and Eta %%%%%%%%%%%%%%%%%%%%
%
% Analyzing the relation between P_FA and eta for the four cases
% with a specified SIR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

sampleSize = 1e8;

numberOfEtaValues = 1000;
maxEta = 4*1e5;
etaValuesCN = [linspace(0.5, 5000, numberOfEtaValues*0.5),linspace(5000, maxEta, numberOfEtaValues*0.5)];
etaValuesCG = [linspace(0.5, 200, numberOfEtaValues*0.2),linspace(200, 1e4, numberOfEtaValues*0.8)];


detectorSigma = 1; % The standard deviation for the detector
clutterSigma = 1; % The standard deviation for the detector
detectorMean = 0;
clutterMean = 0;


SIR = 10; % dB
SIR = 10^(SIR/10);         
alpha = clutterSigma*sqrt(SIR);
theta = 0; 
s = alpha*(cos(theta)+1i*sin(theta)); % signal 

clutterSampleCN = SampleComplexGaussian(sampleSize, clutterMean, clutterSigma); 
clutterSampleCG = SampleCompoundGaussian(sampleSize, clutterMean, clutterSigma); 

sumFA_CNCN = zeros(1, numberOfEtaValues);
sumFA_CNCG = zeros(1, numberOfEtaValues);
sumFA_CGCN = zeros(1, numberOfEtaValues);
sumFA_CGCG = zeros(1, numberOfEtaValues);

% False Alarm (*)
fH1_fa_CNCN = ComplexGaussianPDF(clutterSampleCN, detectorMean + s, detectorSigma); 
fH0_fa_CNCN = ComplexGaussianPDF(clutterSampleCN, detectorMean, detectorSigma);

fH1_fa_CNCG = ComplexGaussianPDF(clutterSampleCG, detectorMean + s, detectorSigma); 
fH0_fa_CNCG = ComplexGaussianPDF(clutterSampleCG, detectorMean, detectorSigma);

fH1_fa_CGCN = CompoundGaussianPDF(clutterSampleCN, detectorMean + s, detectorSigma); 
fH0_fa_CGCN = CompoundGaussianPDF(clutterSampleCN, detectorMean, detectorSigma);

fH1_fa_CGCG = CompoundGaussianPDF(clutterSampleCG, detectorMean + s, detectorSigma); 
fH0_fa_CGCG = CompoundGaussianPDF(clutterSampleCG, detectorMean, detectorSigma);



for iEta=1:numberOfEtaValues
    etaCN = etaValuesCN(iEta);
    etaCG = etaValuesCG(iEta);
    
    sumFA_CNCN(1, iEta) = sum(((fH1_fa_CNCN./fH0_fa_CNCN) > etaCN));

    sumFA_CNCG(1, iEta) = sum(((fH1_fa_CNCG./fH0_fa_CNCG) > etaCN));

    sumFA_CGCN(1, iEta) = sum(((fH1_fa_CGCN./fH0_fa_CGCN) > etaCG));

    sumFA_CGCG(1, iEta) = sum(((fH1_fa_CGCG./fH0_fa_CGCG) > etaCG));
   
    iEta
end

pFalseAlarmCNCN = sumFA_CNCN/sampleSize;
pFalseAlarmCNCG = sumFA_CNCG/sampleSize;
pFalseAlarmCGCN = sumFA_CGCN/sampleSize;
pFalseAlarmCGCG = sumFA_CGCG/sampleSize;



%% Plotting
subplot(1,2,1)
hold on
plot(etaValuesCN,pFalseAlarmCNCN, LineWidth=1.5)
plot(etaValuesCN,pFalseAlarmCNCG, LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('CN-CN', 'CN-CG', location='best')
axis([0, maxEta, 1e-8, 1])

subplot(1,2,2)
set(gca,'ColorOrderIndex',3)
hold on
plot(etaValuesCG,pFalseAlarmCGCN, LineWidth=1.5)
plot(etaValuesCG,pFalseAlarmCGCG, LineWidth=1.5)
set(gca, 'YScale', 'log');
ylabel('P_{FA}'), xlabel('\eta')
legend('CG-CN','CG-CG', location='best')
axis([0, 1e4, 1e-8, 1])




