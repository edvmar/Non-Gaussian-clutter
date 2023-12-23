%%%%%%%%%%%%%% Task2 All known %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear,close all

%% ================== Parameters ========================
% --------- Simulation ---------
sigma = 1;
rMax  = 10*sigma; % kanske större för Kdist? 

numberOfPulses    = 128; % 128
numberOfDistances = 1;  % 100

% --------- Signal ----------- 
omegaD = 0.01;
steeringVector = (exp( 1i*omegaD*(0:numberOfPulses - 1) ))';

SIRs = [0, 3, 5, 7]; % dB 

% ------- Covariance -------- ||| TODO: Seems to be something wrong with Toeplitz. 
epsilon = 1e-10;  % diagonal load
k = 2;
delta   = 1/numberOfPulses^k; % (or 1/numberOfPulses^2)

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta)+ epsilon*eye(numberOfPulses);
L = chol(toeplitzMatrix, 'lower');
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
    F = @(x) 1 - H_nGaussian(abs(x).^2, 0, sigma); % eqn (12)
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

sumFA = zeros(length(SIRs), numberOfEtaValues); 
sumTD = zeros(length(SIRs), numberOfEtaValues);

for iSIR = 1:length(SIRs)
    iSIR
    SIR = 10^(SIRs(iSIR)/10);           
    alpha = sigma*sqrt(SIR);
    signal = alpha*steeringVector;

    % Sampling
    CUTWithoutSignal = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F);
    CUTWithSignal = CUTWithoutSignal + signal; 
    
    % pFA
    q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal,toeplitzMatrixInverse)); 
    q1_H0 = real(MultidimensionalNorm(CUTWithoutSignal-signal,toeplitzMatrixInverse));
    LR_FA = h_n(q1_H0)./h_n(q0_H0);
    
    % pTD
    q0_H1 = real(MultidimensionalNorm(CUTWithSignal,toeplitzMatrixInverse)); 
    q1_H1 = real(MultidimensionalNorm(CUTWithSignal-signal,toeplitzMatrixInverse));
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
%plot([0,1],[0,1])
set(gca, 'XScale', 'log');
xlabel('P_{FA}', FontSize=15), ylabel('P_{TD}',FontSize=15)
legend('SIR = 0', 'SIR = 3', 'SIR = 5', 'SIR = 7', location = 'southeast',FontSize=15)
axis([1e-6, 1, 0, 1])
%savefig('AllKnown_CNK.fig')












