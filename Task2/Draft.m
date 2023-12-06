clear
clc

sampleSize = 100;
rMax  = 5;
sigma = 1;

numberOfPulses    = 10;
numberOfDistances = 8;
epsilon = 1e-6;
delta   = 1/numberOfPulses; % (or 1/numberOfPulses^2)

omegaD  = pi/3; % or 2*pi/numberOfPulsses???
signal = exp( 1i*omegaD*(0:numberOfPulses - 1) )/sqrt(numberOfPulses);

toeplitzMatrix = CalculatePulseCovariance(numberOfPulses, epsilon, delta);
L = chol(toeplitzMatrix + epsilon*eye(numberOfPulses));
toeplitzMatrixInverse = inv(toeplitzMatrix);

% for iSample = 1:sampleSize
rangeBin = SampleComplexGaussianRow(numberOfPulses, rMax, sigma, L);

CPI = zeros(numberOfDistances,numberOfPulses);
for k = 1:numberOfDistances
    rangeBin = SampleComplexGaussianRow(numberOfPulses, rMax, sigma, L);
    CPI(k,:) = rangeBin;
end

CPI;
size(CPI);

% detektor
% P_TD, P_FA

%%
kSignal = 2; % index for range bin with signal in it
CPIwithSignal = CPI;
CPIwithSignal(kSignal,:) = CPIwithSignal(kSignal,:) + signal;


% testing remove later
k0 = 3;
CUT = CPI(k0,:);
CUTsignal = CPI(k0,:) + signal;





% %% Alternativ 1 | Test row with and without signal 
% 
% numberOfEtaValues = 100; 
% etaValues = linspace(0.5,100,numberOfEtaValues);
% 
% sumFA = zeros(1, numberOfEtaValues);
% sumTD = zeros(1, numberOfEtaValues);
% 
% for iEta=1:numberOfEtaValues
%     eta = etaValues(iEta); 
% 
%     for iSample = 1:sampleSize
%     
%         CUT = SampleComplexGaussianRow(numberOfPulses, rMax, sigma, L);
%         CUTsignal = CUT + signal;
%         
%         q0_Clutter = CUT*toeplitzMatrixInverse*CUT';
%         q1_Clutter = (CUT-signal)*toeplitzMatrixInverse*(CUT-signal)';
%         
%         q0_Signal = CUTsignal*toeplitzMatrixInverse*CUTsignal';
%         q1_Signal = (CUTsignal-signal)*toeplitzMatrixInverse*(CUTsignal-signal)';
%         
%         
%         % False Alarm (*)
%         fH1_fa = TailDistributionComplexGaussian(q1_Clutter, numberOfPulses, sigma);          
%         fH0_fa = TailDistributionComplexGaussian(q0_Clutter, numberOfPulses, sigma);          
%         LRT_fa = fH1_fa./fH0_fa;
%         
%         % True Detection (**)
%         fH1_td = TailDistributionComplexGaussian(q1_Signal, numberOfPulses, sigma);           
%         fH0_td = TailDistributionComplexGaussian(q0_Signal, numberOfPulses, sigma);           
%         LRT_td = fH1_td./fH0_td;
% 
%     end
% 
%      % False Alarm (*)
%      sumFA(iEta) = sum((LRT_fa > eta));
%         
%      % True Detection (**)
%      sumTD(iEta) = sum((LRT_td > eta));
% 
% end
% pFalseAlarm = sumFA/sampleSize;
% pDetection  = sumTD/sampleSize;
% 
% % plot(pFalseAlarm, pDetection)














