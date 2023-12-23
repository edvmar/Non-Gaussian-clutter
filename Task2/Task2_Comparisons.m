%%%%%%%%%%%%%% Task2 omega unknown %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

%% ================== Parameters ========================
% --------- Simulation ---------
sampleSize = 1e4; % 2*1e6
sigma = 1;
rMax  = 10*sigma; % kanske större för Kdist? 

numberOfPulses    = 128; % 128
numberOfDistances = 1;  % 100

% --------- Signal ----------- 
SIR = 3; % dB
SIR = 10^(SIR/10);           
alpha = sigma*sqrt(SIR);

% Actual signal
omegaDActual  = 0.01;
steeringVectorActual = (exp( 1i*omegaDActual*(0:numberOfPulses - 1)))';
signalActual = alpha*steeringVectorActual;

% Test signals 
numberOfOmegas = 101;
minOmegaD = 0.0;
maxOmegaD = 0.1;
testOmegaDs  = linspace(minOmegaD, maxOmegaD, numberOfOmegas);

% ------- Covariance ----------------
epsilon = 1e-10;  % diagonal load
k = 2;
delta   = 1/numberOfPulses^k; % (or 1/numberOfPulses^2)

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta)+ epsilon*eye(numberOfPulses);
%toeplitzMatrix = eye(numberOfPulses);
L = chol(toeplitzMatrix,'lower');
toeplitzMatrixInverse = inv(toeplitzMatrix);

% -----  Threshold values ------
numberOfEtaValues = 1000;
etaValues = [linspace(0.01, 1, numberOfEtaValues*0.1),linspace(1, 1e3, numberOfEtaValues*0.75),...
             linspace(2e3, 1e5, numberOfEtaValues*0.1),linspace(1e5, 1e7, numberOfEtaValues*0.05)];


% ------- Distributions ------------
clutterDistribution  = 'CN';  % 'K' or 'CN'
detectorDistribution = 'CN';
nu = 1;


if isequal(clutterDistribution,'CN')
    F = @(x) 1 - H_nGaussian(abs(x).^2, 0, sigma); % eqn (12
elseif isequal(clutterDistribution,'K')
    F = @(x) 1 - H_nKdist(abs(x).^2, 0, sigma, nu); % eqn (12)
else 
    error('Unknown clutter distribution');
end

if isequal(detectorDistribution,'CN')
    h_n = @(x) H_nGaussian(x, numberOfPulses, sigma); % eqn (12
elseif isequal(detectorDistribution,'K')
    h_n = @(x) H_nKdist(x, numberOfPulses, sigma, nu); % eqn (12)
else 
    error('Unknown detector distribution');
end


%% ======================= Simulation ==================================
tic

sumFA = zeros(1, numberOfEtaValues); 
sumTD = zeros(1, numberOfEtaValues);

sumFA_UnknownAlpha = zeros(1, numberOfEtaValues); 
sumTD_UnknownAlpha = zeros(1, numberOfEtaValues);

sumFA_Actual = zeros(1, numberOfEtaValues); 
sumTD_Actual = zeros(1, numberOfEtaValues);

LR_FA = zeros(numberOfOmegas, sampleSize);
LR_TD = zeros(numberOfOmegas, sampleSize);

CUTWithoutSignal = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F);
CUTWithActualSignal = CUTWithoutSignal + signalActual; 


for iOmegaD = 1:numberOfOmegas

    steeringVectorTest = (exp( 1i*testOmegaDs(iOmegaD)*(0:numberOfPulses - 1) )/sqrt(numberOfPulses))';
    steeringVectorNorm = steeringVectorTest'*toeplitzMatrixInverse*steeringVectorTest;
    
    % pFA
    q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal,toeplitzMatrixInverse)); 
    q1_H0 = CalculateQ1WithEstimatedAlpha(CUTWithoutSignal,toeplitzMatrixInverse,...
                                            steeringVectorTest, steeringVectorNorm);
    LR_FA(iOmegaD,:) = h_n(q1_H0)./h_n(q0_H0);
    
    % pTD
    q0_H1 = real(MultidimensionalNorm(CUTWithActualSignal,toeplitzMatrixInverse)); 
    q1_H1 = CalculateQ1WithEstimatedAlpha(CUTWithActualSignal,toeplitzMatrixInverse,...
                                            steeringVectorTest, steeringVectorNorm);
    LR_TD(iOmegaD,:) = h_n(q1_H1)./h_n(q0_H1);
end

% Sort out the max likely omegas 
maxLR_FA = max(LR_FA);
maxLR_TD = max(LR_TD);

% ------- Unknown alpha ---------
steeringVectorNormActual = steeringVectorActual'*toeplitzMatrixInverse*steeringVectorActual;

% pFA
q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal,toeplitzMatrixInverse)); 
q1_H0 = CalculateQ1WithEstimatedAlpha(CUTWithoutSignal,toeplitzMatrixInverse,...
                                            steeringVectorActual, steeringVectorNormActual);
LR_FA_UnknownAlpha = h_n(q1_H0)./h_n(q0_H0);


% pTD
q0_H1 = real(MultidimensionalNorm(CUTWithActualSignal,toeplitzMatrixInverse)); 
q1_H1 = CalculateQ1WithEstimatedAlpha(CUTWithActualSignal,toeplitzMatrixInverse,...
                                            steeringVectorActual, steeringVectorNormActual);
LR_TD_UnknownAlpha = h_n(q1_H1)./h_n(q0_H1);


% ------- All known -----------
% pFA
%q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal,toeplitzMatrixInverse)); % Already calculated
q1_H0 = real(MultidimensionalNorm(CUTWithoutSignal-signalActual,toeplitzMatrixInverse));
LR_FA_Actual = h_n(q1_H0)./h_n(q0_H0);

% pTD
%q0_H1 = real(MultidimensionalNorm(CUTWithActualSignal,toeplitzMatrixInverse)); % Already calculated
q1_H1 = real(MultidimensionalNorm(CUTWithActualSignal-signalActual,toeplitzMatrixInverse));
LR_TD_Actual = h_n(q1_H1)./h_n(q0_H1);


% ----------- Thresholds -----------
for iEta = 1:numberOfEtaValues
    eta = etaValues(iEta);
        
    sumFA(iEta) = sum(maxLR_FA>eta);
    sumTD(iEta) = sum(maxLR_TD>eta);

    sumFA_UnknownAlpha(iEta) = sum((LR_FA_UnknownAlpha>eta));
    sumTD_UnknownAlpha(iEta) = sum((LR_TD_UnknownAlpha>eta));

    sumFA_Actual(iEta) = sum((LR_FA_Actual>eta));
    sumTD_Actual(iEta) = sum((LR_TD_Actual>eta));

end

pFalseAlarm = sumFA/sampleSize;
pDetection  = sumTD/sampleSize;

pFalseAlarmUnknownAlpha = sumFA_UnknownAlpha/sampleSize;
pDetectionUnknownAlpha  = sumTD_UnknownAlpha/sampleSize;

pFalseAlarmActual = sumFA_Actual/sampleSize;
pDetectionActual  = sumTD_Actual/sampleSize;

toc


%% Using stored data 
% 
% load('pActual.mat')
% load('pUnknownAlpha.mat')
% load('pUnknownAlphaOmega.mat')
% load('pUnknownSigma.mat')


%% ============================ Plotting =====================

load('Hej.mat')
figure(1)
hold on
plot(pFalseAlarmActual, pDetectionActual, LineWidth=1.5)
plot(pFalseAlarmUnknownAlpha, pDetectionUnknownAlpha, LineWidth=1.5)
plot(pFalseAlarm, pDetection, LineWidth=1.5)
plot(hejPFA, hejPTD, LineWidth=1.5)
%plot([0,1],[0,1],'k:')
set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend( 'All known', 'Unknown \alpha','Unknown \alpha, \omega','Unknown \alpha, \omega, \Sigma', location = 'best',FontSize=14)
axis([1e-6, 1, 0, 1])











