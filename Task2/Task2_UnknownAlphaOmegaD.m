%%%%%%%%%%%%%% Task2 All known %%%%%%%%%%%%%%
%
% Produces ROC curves for the 1D case where both 
% signal and clutter are known
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

tic 
sampleSize = 1e4;
sigma = 1;
rMax  = 10*sigma; % standardavvikelser kanske större för Kdist? 

numberOfPulses    = 10; % 128
numberOfDistances = 1;  % 100

epsilon = 1e-6;  
k = 0;
delta   = 1/numberOfPulses^k; % (or 1/numberOfPulses^2)
% Seems to be something wrong with Toeplitz. 

numberOfOmegas = 3;
radialVelocity = linspace(1, 100, numberOfOmegas); % m/s
omegaD  = 2*pi*2*radialVelocity/3e8;

SIR = 5; % Loopa flera SIRS sen?
%SIRs = [0, 3, 10, 13]; % dB 
SIR = 10^(SIR/10);           
alpha = sigma*sqrt(SIR);

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta);
%toeplitzMatrix = eye(numberOfPulses);
L = chol(toeplitzMatrix + epsilon*eye(numberOfPulses));
toeplitzMatrixInverse = inv(toeplitzMatrix);
det(toeplitzMatrix)


%%
numberOfEtaValues = 500;
etaValues = linspace(0.001, 10, numberOfEtaValues);

sumFA = zeros(numberOfOmegas, numberOfEtaValues); % Add for other clutters
sumTD = zeros(numberOfOmegas, numberOfEtaValues);

LR_FA = zeros(numberOfOmegas, sampleSize);
LR_TD = zeros(numberOfOmegas, sampleSize);

% Complex Gaussian
F = @(x) 1 - H_nGaussian(abs(x).^2, 0, sigma);  % eqn (12)    % Clutter dist
h_n = @(x) H_nGaussian(x, numberOfPulses, sigma);             % Detector dist

% complex K distribution
nu = 0.01;
%F = @(x) 1 - H_nKdist(abs(x).^2, 0, sigma, nu);  % eqn (12)  % Clutter dist
%h_n = @(x) H_nKdist(x, numberOfPulses, sigma, nu);           % Detector dist




% Sampling.. gör snabbare senare 
for iSample = 1:sampleSize

    CUTWithoutTheSignal = Sampling(numberOfPulses, rMax, sigma, L, F);
    for iOmegaD = 1:numberOfOmegas
        steeringVector = exp( 1i*omegaD(iOmegaD)*(0:numberOfPulses - 1) )/sqrt(numberOfPulses);
        steeringVectorNorm = steeringVector*toeplitzMatrixInverse*steeringVector';
        signal = alpha*steeringVector; 
        CUTWithSignal = CUTWithoutTheSignal + signal; 
        
        % pFA
        q0_H0 = real(CUTWithoutTheSignal*toeplitzMatrixInverse*CUTWithoutTheSignal');
        q1_H0 = CalculateQ1WithEstimatedAlpha(CUTWithoutTheSignal,toeplitzMatrixInverse,...
                                                steeringVector, steeringVectorNorm);
        LR_FA(iOmegaD, iSample) = h_n(q1_H0)/h_n(q0_H0);
        
    
        % pTD
        q0_H1 = real(CUTWithSignal*toeplitzMatrixInverse*CUTWithSignal');
        q1_H1 = CalculateQ1WithEstimatedAlpha(CUTWithSignal,toeplitzMatrixInverse,...
                                                steeringVector, steeringVectorNorm);
        LR_TD(iOmegaD, iSample) = h_n(q1_H1)/h_n(q0_H1);


    end
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

%%
hold on
for iOmegaD = 1:numberOfOmegas
   plot(pFalseAlarm(iOmegaD,:), pDetection(iOmegaD, :), LineWidth=1.5)
end
plot([0,1],[0,1])
%set(gca, 'XScale', 'log');
xlabel('P_{FA}'), ylabel('P_{TD}')
%legend('SIR = 0', 'SIR = 3', 'SIR = 10', 'SIR = 13', location = 'southeast')
%axis([1e-7, 1, 0, 1])











