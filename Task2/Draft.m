
sampleSize = 100;
rMax = 5;
sigma = 1;

numberOfPulses = 10;
numberOfDistances = 8;
epsilon = 1e-6;
delta = 1/numberOfPulses; %(or 1/numberOfPulses^2)

L = CalculatePulseCovariance(numberOfPulses, epsilon, delta);

% for iSample = 1:sampleSize
clutterRow = SampleComplexGaussianRow(numberOfPulses, rMax,sigma,L);

CPI = zeros(numberOfDistances,numberOfPulses);
for k = 1:numberOfDistances
    clutterRow = SampleComplexGaussianRow(numberOfPulses, rMax,sigma,L);
    CPI(k,:) = clutterRow;
end

CPI
size(CPI)

% detektor
% P_TD, P_FA

