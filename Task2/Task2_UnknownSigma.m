%%%%%%%%%%%%%% Task2 Sigma Unknown %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

%% ================== Parameters ========================
% --------- Simulation ---------
sampleSize = 1e4;
sigma = 1;
rMax  = 10*sigma; % kanske större för Kdist? 

numberOfPulses    = 128; % 128
numberOfDistances = 100;  % 100

% --------- Signal ----------- 
SIRs = [0, 3, 5, 7]; % dB 

% Actual signal
omegaDActual  = 0.01;
steeringVectorActual = (exp( 1i*omegaDActual*(0:numberOfPulses - 1)))';

% Test signals 
numberOfOmegas = 101;
minOmegaD = 0.0;
maxOmegaD = 0.1;
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
etaValues = [linspace(0.1, 500, numberOfEtaValues*0.3),linspace(500, 1000000, numberOfEtaValues*0.7)];


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

sumFA = zeros(length(SIRs), numberOfEtaValues); 
sumTD = zeros(length(SIRs), numberOfEtaValues);
signalRow = numberOfDistances-1;

for iSIR = 1:length(SIRs)
    iSIR
    
    SIR = 10^(SIRs(iSIR)/10);           
    alpha = sigma*sqrt(SIR);
    signalActual = alpha*steeringVectorActual;


    for iSample=1:sampleSize
        % Sampling
        CPI = Sampling(numberOfPulses, numberOfDistances, rMax, sigma, L, F)';
        
        CUTWithoutSignal = CPI(signalRow,:)';
        CUTWithActualSignal = CUTWithoutSignal + signalActual;
        CPIWithoutSignal = CPI;
        CPIWithoutSignal(signalRow,:) = [];
        
        covarianceEstimate = 1/(numberOfDistances-1)*(CPIWithoutSignal')*CPIWithoutSignal; 
        covarianceEstimate = real(covarianceEstimate);
        covarianceEstimateInverse = inv(covarianceEstimate + epsilon*eye(numberOfPulses));
        % Include CUT in covariance estimate? 
        
        %covarianceEstimateInverse = toeplitzMatrixInverse;

        LR_FA = zeros(1,numberOfOmegas);
        LR_TD = zeros(1,numberOfOmegas);

        for iOmegaD = 1:numberOfOmegas
        
            steeringVectorTest = (exp( 1i*testOmegaDs(iOmegaD)*(0:numberOfPulses - 1) )/sqrt(numberOfPulses))';
            steeringVectorNorm = steeringVectorTest'*covarianceEstimateInverse*steeringVectorTest;
            
            % pFA
            q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal,covarianceEstimateInverse)); 
            q1_H0 = CalculateQ1WithEstimatedAlpha(CUTWithoutSignal,covarianceEstimateInverse,...
                                                    steeringVectorTest, steeringVectorNorm);
            LR_FA(iOmegaD) = h_n(q1_H0)./h_n(q0_H0);
            
            % pTD
            q0_H1 = real(MultidimensionalNorm(CUTWithActualSignal,covarianceEstimateInverse)); 
            q1_H1 = CalculateQ1WithEstimatedAlpha(CUTWithActualSignal,covarianceEstimateInverse,...
                                                    steeringVectorTest, steeringVectorNorm);
            LR_TD(iOmegaD) = h_n(q1_H1)./h_n(q0_H1);
        end
    
        % Sort out the max likely omegas 
        maxLR_FA = max(LR_FA);
        maxLR_TD = max(LR_TD);
        
        sumFA(iSIR,:) = sumFA(iSIR,:) + (maxLR_FA>etaValues);
        sumTD(iSIR,:) = sumTD(iSIR,:) + (maxLR_TD>etaValues);
   
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
xlabel('P_{FA}', FontSize=15), ylabel('P_{TD}',FontSize=15)
legend('SIR = 0', 'SIR = 3', 'SIR = 5', 'SIR = 7', location = 'southeast',FontSize=15)
axis([1e-7, 1, 0, 1])











