%%%%%%%%%%%%%%%%%%% Testing Distributions %%%%%%%%%%%%%%%%%%%%
%
% Produces figures to test the feasibility of the 
% implemented distributions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% Testing the sample
sampleSize = 1e3;
sigma = 1;
sample1 = SampleComplexGaussian(sampleSize, 0, sigma);
sample2 = SampleCompoundGaussian(sampleSize, 0, sigma);

figure(1)
subplot(1,2,1)
histogram(abs(sample1))
subplot(1,2,2)
histogram(abs(sample2))

figure(2)
subplot(1,2,1)
histogram2(real(sample1), imag(sample1))
subplot(1,2,2)
histogram2(real(sample2), imag(sample2))

figure(3)
subplot(1,2,1)
plot(real(sample1), imag(sample1), 'bo')
axis equal
axis([-4 4 -4 4])
subplot(1,2,2)
plot(real(sample2), imag(sample2), 'ro')
axis equal
axis([-4 4 -4 4])

figure(4)
[X,Y] = meshgrid(-5:0.1:5,-5:0.1:5);
subplot(1,2,1)
Z1 = ComplexGaussianPDF(X+1i*Y, 0, sigma)*sampleSize;
surf(X,Y,Z1)
subplot(1,2,2)
Z2 = CompoundGaussianPDF(X+1i*Y, 0, sigma)*sampleSize;
surf(X,Y,Z2)



%% Testing the PDF
clc,clear
figure(5)
sigma = 1;
x = linspace(-5,5);
y1 = ComplexGaussianPDF(x, 0, sigma);
y2 = CompoundGaussianPDF(x, 0, sigma);
hold on
plot(x, y1, 'b', LineWidth=1.5)
plot(x, y2, 'r', LineWidth=1.5)
legend('CN', 'CG')

% f=@(x) exp(-sqrt(2)*abs(x)/sigma)./(sqrt(2)*pi*abs(x));
% plot(x, f(x),'g',Linewidth = 1.5)

% fun=@(x,y) exp(-sqrt(2)*sqrt(x.^2+y.^2)/sigma)./(sqrt(2)*pi*(x.^2+y.^2));
% integral2(fun ,-1000, 1000, -1000, 1000);


%% Calculate variances
clear
clc

sigma = 3;
fComp = @(a,b)besselk(0, 2*sqrt(a.^2+b.^2)./sigma)/(2*pi*(sigma/2)^2);
fVar = @(a,b) (a.^2+b.^2).*besselk(0, 2*sqrt(a.^2+b.^2)./sigma)/(2*pi*(sigma/2)^2);
integral2(fComp ,-100, 100, -100, 100);
Var = integral2(fVar, -100, 100, -100, 100)
sigmaNy = sqrt(Var);


% sigmaCN = 2*sigma = sqrt(variance)

%
f1dim = @(x) besselk(0, abs(x)/sigma)/(pi*sigma);
f1dimVar = @(x) x.^2.*besselk(0, abs(x)/sigma)/(pi*sigma);
integral(f1dim,-1000,1000);
Var1dim = integral(f1dimVar,-1000,1000);
sigma1dim = sqrt(Var1dim);


sampleSize = 1e5;
sigmaPrior = sigma;
sample = SampleCompoundGaussian(sampleSize, 0, sigmaPrior);

var(real(sample))+var(imag(sample))











