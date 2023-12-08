%%%%%%%%%%%%%% Task2 UnkownAlpha %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known unkown alpha
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

% ================== Parameters ========================
% --------- Simulation ---------
sampleSize = 1e4;
sigma = 1;
rMax  = 10*sigma; % kanske större för Kdist? 

numberOfPulses    = 10; % 128
numberOfDistances = 1;  % 100

% --------- Signal ----------- 
radialVelocity = 100; % m/s
omegaD  = 2*pi*2*radialVelocity/3e8; % Doppler Freq
steeringVector = (exp( 1i*omegaD*(0:numberOfPulses - 1) )/sqrt(numberOfPulses))';

SIR = 5; % Loopa flera SIRS sen?
%SIRs = [0, 3, 10, 13]; % dB 
SIR = 10^(SIR/10);           
alpha = sigma*sqrt(SIR);
signal = alpha*steeringVector;

% ------- Covariance -------- ||| TODO: Seems to be something wrong with Toeplitz. 
epsilon = 1e-6;  % diagonal load
k = 0;
delta   = 1/numberOfPulses^k; % (or 1/numberOfPulses^2)

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta);
%toeplitzMatrix = eye(numberOfPulses);
L = chol(toeplitzMatrix + epsilon*eye(numberOfPulses));
toeplitzMatrixInverse = inv(toeplitzMatrix);
det(toeplitzMatrix)

% -----  Threshold values ------
numberOfEtaValues = 500;
etaValues = linspace(0.001, 100, numberOfEtaValues);

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

steeringVectorNorm = steeringVector'*toeplitzMatrixInverse*steeringVector;

%% ======================= Simulation ==================================
tic
sumFA = zeros(1, numberOfEtaValues); % Add for other clutters
sumTD = zeros(1, numberOfEtaValues);

CUTWithoutSignal = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F);
CUTWithSignal = CUTWithoutSignal + signal; 

% pFA
q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal,toeplitzMatrixInverse)); 
q1_H0 = CalculateQ1WithEstimatedAlpha(CUTWithoutSignal,toeplitzMatrixInverse,...
                                            steeringVector, steeringVectorNorm);
LR_FA = h_n(q1_H0)./h_n(q0_H0);


% pTD
q0_H1 = real(MultidimensionalNorm(CUTWithSignal,toeplitzMatrixInverse)); 
q1_H1 = CalculateQ1WithEstimatedAlpha(CUTWithSignal,toeplitzMatrixInverse,...
                                            steeringVector, steeringVectorNorm);
LR_TD = h_n(q1_H1)./h_n(q0_H1);

for iEta = 1:numberOfEtaValues

        eta = etaValues(iEta);

        sumFA(1, iEta) = sum((LR_FA>eta));
        sumTD(1, iEta) = sum((LR_TD>eta));
end

pFalseAlarm = sumFA/sampleSize;
pDetection = sumTD/sampleSize;

toc

%% ============================ Plotting =====================
hold on
%for iSIR = 1:length(SIRs)
%    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
%end
plot(pFalseAlarm, pDetection, LineWidth = 1.5)
plot([0,1],[0,1])
%set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
%legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location = 'southeast')
%axis([1e-7, 1, 0, 1])











