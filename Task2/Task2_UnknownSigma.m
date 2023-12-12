%%%%%%%%%%%%%% Task2 All known %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

%% ================== Parameters ========================
% --------- Simulation ---------
sampleSize = 1000;
sigma = 1;
rMax  = 10*sigma; % kanske större för Kdist? 

numberOfPulses    = 128; % 128
numberOfDistances = 100;  % 100

% --------- Signal ----------- 
radialVelocity = 100; % m/s
omegaD  = 2*pi*2*radialVelocity/3e8; % Doppler Freq
steeringVector = (exp( 1i*omegaD*(0:numberOfPulses - 1) )/sqrt(numberOfPulses))';

SIR = 10; % Loopa flera SIRS sen?
%SIRs = [0, 3, 10, 13]; % dB 
SIR = 10^(SIR/10);           
alpha = sigma*sqrt(SIR);
signal = alpha*steeringVector;

% ------- Covariance -------- ||| TODO: Seems to be something wrong with Toeplitz. 
epsilon = 1e-10;  % diagonal load
k = 0;
delta   = 1/numberOfPulses^k; % (or 1/numberOfPulses^2)

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta);
%toeplitzMatrix = eye(numberOfPulses);
L = chol(toeplitzMatrix + epsilon*eye(numberOfPulses));
toeplitzMatrixInverse = inv(toeplitzMatrix);

% -----  Threshold values ------
numberOfEtaValues = 500;
etaValues = [linspace(1,100, numberOfEtaValues*0.1),linspace(100, 100000, numberOfEtaValues*0.9)];

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
sumFA = zeros(1, numberOfEtaValues); % Add for other clutters
sumTD = zeros(1, numberOfEtaValues);
row_where_signal = numberOfDistances-1;

for i=1:sampleSize
% Sampling
    CPI = Sampling(numberOfPulses, numberOfDistances, rMax, sigma, L, F)';
    CUTWithoutSignal = CPI(row_where_signal,:);
    CUTWithSignal = CUTWithoutSignal + signal';
    CPIWithoutSignal = CPI;
    CPIWithoutSignal(row_where_signal,:) = [];
    
    covarianceEst = 1/(numberOfDistances-1)*(CPIWithoutSignal')*CPIWithoutSignal;
    covarianceEst = real(covarianceEst);
    covarianceEstInverse = inv(covarianceEst);
    
    steeringVectorNorm = steeringVector'*covarianceEstInverse*steeringVector;
    
    % pFA
    q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal',covarianceEstInverse)); 
    q1_H0 = CalculateQ1WithEstimatedAlpha(CUTWithoutSignal',covarianceEstInverse,...
                                            steeringVector, steeringVectorNorm);
    LR_FA = h_n(q1_H0)./h_n(q0_H0);
    
    % pTD
    q0_H1 = real(MultidimensionalNorm(CUTWithSignal',covarianceEstInverse)); 
    q1_H1 = CalculateQ1WithEstimatedAlpha(CUTWithSignal',covarianceEstInverse,...
                                            steeringVector, steeringVectorNorm);
    
    LR_TD = h_n(q1_H1)./h_n(q0_H1);

    for iEta = 1:numberOfEtaValues
    
        eta = etaValues(iEta);
    
        sumFA(iEta) = sumFA(iEta) + sum((LR_FA>eta));
        sumTD(iEta) = sumTD(iEta) + sum((LR_TD>eta));
    end
end
pFalseAlarm = sumFA/sampleSize;
pDetection = sumTD/sampleSize;

toc

%% ============================ Plotting =====================
hold on
%for iSIR = 1:length(SIRs)
%    plot(pFalseAlarm(iSIR,:), pDetection(iSIR, :), LineWidth=1.5)
%end
plot(pFalseAlarm, pDetection, ':',LineWidth = 1.5)
plot([0,1],[0,1],'k:')
%set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend('SIR = 10', location = 'best')
%legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location = 'southeast')
%axis([1e-7, 1, 0, 1])











