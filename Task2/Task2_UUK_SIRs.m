%%%%%%%%%%%%%% Task2 omega unknown %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear%, close all

%% ================== Parameters ========================
% --------- Simulation ---------
sampleSize = 1e4;
sigma = 1;
rMax  = 10*sigma; % kanske större för Kdist? 

numberOfPulses    = 128; % 128
numberOfDistances = 1;  % 100

% --------- Signal ----------- 
SIRs = [0, 1, 3, 5]; % dB 

% Actual signal
omegaDActual  = 0.01;
steeringVectorActual = (exp( 1i*omegaDActual*(0:numberOfPulses - 1)))';

% Test signals 
numberOfOmegas = 100;
minOmegaD = 0.005;
maxOmegaD = 0.05;
testOmegaDs  = linspace(minOmegaD, maxOmegaD, numberOfOmegas);

% ------- Covariance ----------------
epsilon = 1e-10;  % diagonal load
k = 2;
delta   = 0.5/numberOfPulses^k; % (or 1/numberOfPulses^2)

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta)+ epsilon*eye(numberOfPulses);
%toeplitzMatrix = eye(numberOfPulses);
L = chol(toeplitzMatrix,'lower');
toeplitzMatrixInverse = inv(toeplitzMatrix);

% -----  Threshold values ------
numberOfEtaValues = 1000;
etaValues = [linspace(0.1, 100, numberOfEtaValues*0.3),linspace(100, 10000, numberOfEtaValues*0.7)];


% ------- Distributions ------------
clutterDistribution  = 'K';  % 'K' or 'CN'
detectorDistribution = 'K';
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

sumFA = zeros(length(SIRs), numberOfEtaValues); 
sumTD = zeros(length(SIRs), numberOfEtaValues);

for iSIR = 1:length(SIRs)
    iSIR
    SIR = 10^(SIRs(iSIR)/10);           
    alpha = sigma*sqrt(SIR);
    signalActual = alpha*steeringVectorActual;

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

    % ----------- Thresholds -----------
    for iEta = 1:numberOfEtaValues
        eta = etaValues(iEta);
            
        sumFA(iSIR,iEta) = sum(maxLR_FA>eta);
        sumTD(iSIR,iEta) = sum(maxLR_TD>eta);
    
    end
end

pFalseAlarm = sumFA/sampleSize;
pDetection  = sumTD/sampleSize;


toc

%% ============================ Plotting =====================
hold on
for iSIR = 1:length(SIRs)
    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
end
%plot([0,1],[0,1])
set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend('SIR = 0', 'SIR = 1', 'SIR = 3', 'SIR = 5', location = 'best')
axis([1e-7, 1, 0, 1])










