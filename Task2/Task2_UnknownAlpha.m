%%%%%%%%%%%%%% Task2 UnkownAlpha %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where clutter is known
% but the signal strength alpha is unknown
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

% ================== Parameters ========================
% --------- Simulation ---------
sampleSize = 2*1e6;
sigma = 1;
rMax  = 10*sigma; 

numberOfPulses    = 128; 
numberOfDistances = 1;  % Only need one CUT

% --------- Signal ----------- 
omegaD = 0.01;
steeringVector = (exp( 1i*omegaD*(0:numberOfPulses - 1) ))';

SIRs = [0, 3, 5, 7]; % dB 


% ------- Covariance -------- 
epsilon = 1e-10;  % diagonal load
delta   = 1/numberOfPulses^2;

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta)+ epsilon*eye(numberOfPulses);
L = chol(toeplitzMatrix, 'lower');
toeplitzMatrixInverse = inv(toeplitzMatrix);

% -----  Threshold values ------
numberOfEtaValues = 1000;
etaValues = [linspace(1, 100, numberOfEtaValues*0.1),linspace(110, 1e4, numberOfEtaValues*0.5),...
             linspace(2e4, 1e7, numberOfEtaValues*0.2),linspace(1e7, 1e10, numberOfEtaValues*0.2)];


% ------- Distributions ------------
clutterDistribution  = 'CN';  % 'K' or 'CN'
detectorDistribution = 'CN';
nu = 1;


if isequal(clutterDistribution,'CN')
    F = @(x) 1 - H_nGaussian(abs(x).^2, 0, sigma); % CDF
elseif isequal(clutterDistribution,'K')
    F = @(x) 1 - H_nKdist(abs(x).^2, 0, sigma, nu); 
else 
    error('Unknown clutter distribution');
end

if isequal(detectorDistribution,'CN')
    h_n = @(x) H_nGaussian(x, numberOfPulses, sigma); % TailDistribution
elseif isequal(detectorDistribution,'K')
    h_n = @(x) H_nKdist(x, numberOfPulses, sigma, nu);
else 
    error('Unknown detector distribution');
end

%% ======================= Simulation ==================================
tic

steeringVectorNorm = steeringVector'*toeplitzMatrixInverse*steeringVector;

sumFA = zeros(length(SIRs), numberOfEtaValues);
sumTD = zeros(length(SIRs), numberOfEtaValues);

for iSIR = 1:length(SIRs)
    iSIR
    SIR = 10^(SIRs(iSIR)/10);           
    alpha = sigma*sqrt(SIR);
    signal = alpha*steeringVector;

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
    
            sumFA(iSIR, iEta) = sum((LR_FA>eta));
            sumTD(iSIR, iEta) = sum((LR_TD>eta));
    end
end

pFalseAlarm = sumFA/sampleSize;
pDetection = sumTD/sampleSize;

toc

%% ============================ Plotting =====================
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
end
set(gca, 'XScale', 'log');
xlabel('P_{FA}', FontSize=15), ylabel('P_{TD}',FontSize=15)
legend('SIR = 0', 'SIR = 3', 'SIR = 5', 'SIR = 7', location = 'southeast',FontSize=15)
axis([1e-6, 1, 0, 1])











