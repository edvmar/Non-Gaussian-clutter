%%%%%%%%%%%%%% Task2 omega unknown %%%%%%%%%%%%%%
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
numberOfDistances = 1;  % 100

% --------- Signal ----------- 
SIR = 20; 
SIR = 10^(SIR/10);           
alpha = sigma*sqrt(SIR);

% Actual signal
actualRadialVelocity = 25; %m/s 
omegaD  = 2*pi*2*actualRadialVelocity/3e8;
steeringVector = (exp( 1i*omegaD*(0:numberOfPulses - 1) )/sqrt(numberOfPulses))';
signal = alpha*steeringVector;

% Test signals 
numberOfOmegas = 5;
maxRadialVelocity = 100; % m/s TODO: borde vara 100 m/s
radialVelocities = [linspace(1, maxRadialVelocity, numberOfOmegas-1)]; % m/s
radialVelocities = [radialVelocities, actualRadialVelocity];
omegaDs  = 2*pi*2*radialVelocities/3e8;


% ------- Covariance -------- ||| TODO: Seems to be something wrong with Toeplitz. 
epsilon = 1e-10;  % diagonal load
k = 2;
delta   = 1/numberOfPulses^k; % (or 1/numberOfPulses^2)

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta)+ epsilon*eye(numberOfPulses);
%toeplitzMatrix = eye(numberOfPulses);
L = chol(toeplitzMatrix,'lower');
toeplitzMatrixInverse = inv(toeplitzMatrix);

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




%% ======================= Simulation ==================================
tic

sumFA = zeros(numberOfOmegas, numberOfEtaValues); % Add for other clutters
sumTD = zeros(numberOfOmegas, numberOfEtaValues);

LR_FA = zeros(numberOfOmegas, sampleSize);
LR_TD = zeros(numberOfOmegas, sampleSize);

CUTWithoutSignal = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F);
CUTWithActualSignal = CUTWithoutSignal + signal; 


for iOmegaD = 1:numberOfOmegas
    steeringVector = (exp( 1i*omegaDs(iOmegaD)*(0:numberOfPulses - 1) )/sqrt(numberOfPulses))';
    steeringVectorNorm = steeringVector'*toeplitzMatrixInverse*steeringVector;
    
    % pFA
    q0_H0 = real(MultidimensionalNorm(CUTWithoutSignal,toeplitzMatrixInverse)); 
    q1_H0 = CalculateQ1WithEstimatedAlpha(CUTWithoutSignal,toeplitzMatrixInverse,...
                                            steeringVector, steeringVectorNorm);
    LR_FA(iOmegaD,:) = h_n(q1_H0)./h_n(q0_H0);
    
    % pTD
    q0_H1 = real(MultidimensionalNorm(CUTWithActualSignal,toeplitzMatrixInverse)); 
    q1_H1 = CalculateQ1WithEstimatedAlpha(CUTWithActualSignal,toeplitzMatrixInverse,...
                                            steeringVector, steeringVectorNorm);
    LR_TD(iOmegaD,:) = h_n(q1_H1)./h_n(q0_H1);
end


for iEta = 1:numberOfEtaValues
    eta = etaValues(iEta);
    for iOmegaD = 1:numberOfOmegas

        sumFA(iOmegaD, iEta) = sum((LR_FA(iOmegaD,:)>eta));
        sumTD(iOmegaD, iEta) = sum((LR_TD(iOmegaD,:)>eta));

    end
end

pFalseAlarm = sumFA/sampleSize;
pDetection  = sumTD/sampleSize;

toc

%% ============================ Plotting =====================
hold on
for iOmegaD = 1:numberOfOmegas-1
   plot(pFalseAlarm(iOmegaD,:), pDetection(iOmegaD, :), LineWidth=1.5)
end
plot(pFalseAlarm(numberOfOmegas,:), pDetection(numberOfOmegas, :), 'k--', LineWidth=1.5) % NOTE: Last index is the actual
plot([0,1],[0,1])
%set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
legend('w1', 'w2', 'w3', 'w4','wActual', location = 'southeast')
%axis([1e-7, 1, 0, 1])











