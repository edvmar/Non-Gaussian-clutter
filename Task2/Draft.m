clear
clc

sampleSize = 100;
rMax  = 5;
sigma = 1;

numberOfPulses    = 10;
numberOfDistances = 8;
epsilon = 1e-6;
delta   = 1/numberOfPulses; % (or 1/numberOfPulses^2)

omegaD  = pi/3; % or 2*pi/numberOfPulses???
signal = exp( 1i*omegaD*(1:numberOfPulses - 1) )/sqrt(numberOfPulses);

L = CalculatePulseCovariance(numberOfPulses, epsilon, delta);

% for iSample = 1:sampleSize
rangeBin = SampleComplexGaussianRow(numberOfPulses, rMax,sigma,L);

CPI = zeros(numberOfDistances,numberOfPulses);
for k = 1:numberOfDistances
    rangeBin = SampleComplexGaussianRow(numberOfPulses, rMax,sigma,L);
    CPI(k,:) = rangeBin;
end

kSignal = 2; % index for range bin with signal in it
CPIwithSignal = CPI;
CPIwithSignal(kSignal,:) = CPIwithSignal(kSignal,:) + signal;

CPI;
size(CPI);

% detektor
% P_TD, P_FA

