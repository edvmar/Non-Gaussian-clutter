clear
clc


%% ================== Parameters ========================
% --------- Simulation ---------
sampleSize = 4;
sigma = 1;
rMax  = 3.3*sigma; % kanske större för Kdist? 

numberOfPulses    = 6; % 128
numberOfDistances = 5;  % 100

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

% Table for cdf x,y values
domain = linspace(0, rMax, 7); % Byt till inget rMax beroende?
range = F(domain)';
rangeMatrix = repmat(range, 1, numberOfPulses, numberOfDistances, sampleSize);
rangeMatrix = permute(rangeMatrix, [2 3 4 1]); % to get correct dimensions to calc rangeDifference

uniformSample = rand(numberOfPulses, numberOfDistances, sampleSize);

rangeDifference = rangeMatrix-uniformSample;
rangeDifference = permute(rangeDifference, [4, 1, 2, 3]); % to get correct dimensions for the min operation

[~, index] = min(abs(rangeDifference));
radius = domain(index);
radius = squeeze(radius); % removes the redundant first dimension of radius

thetas = rand(1,numberOfDistances)*2*pi;
xSample = exp(thetas*1i).*radius; 

CPIsamples_test = pagemtimes(L,xSample); % pagewise matrix multiplication
                                         % each page corresponds to a CPI
                                         % matrix sample
