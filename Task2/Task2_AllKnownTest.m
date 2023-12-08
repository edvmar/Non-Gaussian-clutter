%%%%%%%%%%%%%% Task2 All known %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

sampleSize = 10;

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

toeplitzMatrix = CalculatePulseCovariance(numberOfPulses, delta);
L = chol(toeplitzMatrix + epsilon*eye(numberOfPulses));
toeplitzMatrixInverse = inv(toeplitzMatrix);


numberOfEtaValues = 5;
etaValues = linspace(0.001, 100, numberOfEtaValues);

sumFA = zeros(1, numberOfEtaValues); % Add for other clutters
sumTD = zeros(1, numberOfEtaValues);



q0_H0 = zeros(1, sampleSize);
q1_H0 = zeros(1, sampleSize);

% Sampling.. gör snabbare senare 
for i = 1:sampleSize
    
    CUTnoSignal = SampleComplexGaussianRow(numberOfPulses, rMax, sigma, L);
    CUTsignal = CUTnoSignal + signal;
    
    % pFA
    q0_H0(i) = real(CUTnoSignal*toeplitzMatrixInverse*CUTnoSignal');
    q1_H0(i) = real((CUTnoSignal-signal)*toeplitzMatrixInverse*(CUTnoSignal-signal)');
    
    
    % pTD
    q0_H1 = real(CUTsignal*toeplitzMatrixInverse*CUTsignal');
    q1_H1 = real((CUTsignal-signal)*toeplitzMatrixInverse*(CUTsignal-signal)');
    LR_TD(i) = TailDistributionComplexGaussian(q1_H1, numberOfPulses, sigma)/...
                    TailDistributionComplexGaussian(q0_H1, numberOfPulses, sigma);

end

LR_FA(i) = TailDistributionComplexGaussian(q1_H0, numberOfPulses, sigma)./...
                   TailDistributionComplexGaussian(q0_H0, numberOfPulses, sigma);


for iEta = 1:numberOfEtaValues

        eta = etaValues(iEta);

        sumFA(1, iEta) = sum((LR_FA>eta));
        sumTD(1, iEta) = sum((LR_TD>eta));

end

pFalseAlarm = sumFA/sampleSize;
pDetection = sumTD/sampleSize;

%%
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










