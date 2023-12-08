%%%%%%%%%%%%% Sample Testing %%%%%%%%%%%%%%%%%
%
% Trying to implement the sampling with functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
sampleSize = 1e4;
sigma = 1;
nu = 1;
rMax = 10*sigma;
numberOfPulses = 1;
epsilon = 1e-8;
delta = 1/numberOfPulses;

toeplitzMatrix = CalculateToeplitzMatrix(numberOfPulses, delta);
L = chol(toeplitzMatrix + epsilon*eye(numberOfPulses));


F_CG = @(x) 1 - H_nGaussian(abs(x).^2, 0, sigma);  % eqn (12)    % Clutter dist
h_n_CG = @(x) H_nGaussian(x, numberOfPulses, sigma);             % Detector dist

% complex K distribution
F_K = @(x) 1 - H_nKdist(abs(x).^2, 0, sigma, nu);  % eqn (12)  % Clutter dist
h_n_K = @(x) H_nKdist(x, numberOfPulses, sigma, nu); 


sample1 = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F_CG);
sample2 = Sampling(numberOfPulses, sampleSize, rMax, sigma, L, F_K);

%%
figure(1)
subplot(1,2,1)
histogram2(real(sample1), imag(sample1))
subplot(1,2,2)
histogram2(real(sample2), imag(sample2))

figure(2)
subplot(1,2,1)
plot(real(sample1), imag(sample1), 'bo')
axis equal
%rMax=2*rMax;
axis([-rMax rMax -rMax rMax])
subplot(1,2,2)
plot(real(sample2), imag(sample2), 'ro')
axis equal
axis([-rMax rMax -rMax rMax])
%
x = linspace(-rMax,rMax);
sigma1D = sigma/sqrt(2);
f = @(x) 1/(sqrt(2*pi*sigma1D^2))*exp(-x.^2/(2*sigma1D^2));
figure(3)

subplot(1,2,1)
hold on
histogram(real(sample1), 'Normalization','pdf')
plot(x,f(x))

subplot(1,2,2)
hold on
histogram(real(sample2), 'Normalization','pdf')
plot(x,f(x))
