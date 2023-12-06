%%%%%%%%%%%%%% Task2 All known %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all


sampleSize = 1e4;
sigma = 1;
rMax  = 10*sigma; % standardavvikelser

numberOfPulses    = 10; % 128
numberOfDistances = 8;  % 100

epsilon = 1e-10;      
delta   = 1/numberOfPulses; % (or 1/numberOfPulses^2)

omegaD  = 1e-7; % Check later.. 
steeringVector = exp( 1i*omegaD*(0:numberOfPulses - 1) )/sqrt(numberOfPulses);

SIR = 10; % Loopa flera SIRS sen?
%SIRs = [0, 3, 10, 13]; % dB 
SIR = 10^(SIR/10);           
alpha = sigma*sqrt(SIR);
signal = alpha*steeringVector;

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta);
L = chol(toeplitzMatrix + epsilon*eye(numberOfPulses));
toeplitzMatrixInverse = inv(toeplitzMatrix);


numberOfEtaValues = 5000;
etaValues = linspace(0.001, 100, numberOfEtaValues);

sumFA = zeros(1, numberOfEtaValues); % Add for other clutters
sumTD = zeros(1, numberOfEtaValues);

LR_FA = zeros(1, sampleSize);
LR_TD = zeros(1, sampleSize);

% Complex Gaussian
F = @(x) 1 - H_nGaussian(abs(x).^2, 0, sigma);  % eqn (12) 
h_n = @(x) H_nGaussian(x, numberOfPulses, sigma);

% complex K distribution

% Sampling.. gör snabbare senare 
for i = 1:sampleSize
    
    CUTWithoutTheSignal = Sampling(numberOfPulses, rMax, sigma, L, F); % Feels better having two different rows
    CUTWithSignal = Sampling(numberOfPulses, rMax, sigma, L, F) + signal; 
    
    % pFA
    q0_H0 = real(CUTWithoutTheSignal*toeplitzMatrixInverse*CUTWithoutTheSignal');
    q1_H0 = real((CUTWithoutTheSignal-signal)*toeplitzMatrixInverse*(CUTWithoutTheSignal-signal)');
    LR_FA(i) = h_n(q1_H0)/h_n(q0_H0);
    

    % pTD
    q0_H1 = real(CUTWithSignal*toeplitzMatrixInverse*CUTWithSignal');
    q1_H1 = real((CUTWithSignal-signal)*toeplitzMatrixInverse*(CUTWithSignal-signal)');
    LR_TD(i) = h_n(q1_H1)/h_n(q0_H1);

end


for iEta = 1:numberOfEtaValues

        eta = etaValues(iEta);

        sumFA(1, iEta) = sum((LR_FA>eta));
        sumTD(1, iEta) = sum((LR_TD>eta));
end

pFalseAlarm = sumFA/sampleSize;
pDetection = sumTD/sampleSize;

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











