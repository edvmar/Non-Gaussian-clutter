clear
clc

%% Testing the sample

sample1 = SampleComplexGaussian(10^3, 0, 1);
sample2 = SampleCompoundGaussian(10^3, 0, 1);

figure(1)
subplot(1,2,1)
histogram(abs(sample1))
subplot(1,2,2)
histogram(abs(sample2))


figure(2)
subplot(1,2,1)
plot(real(sample1), imag(sample1), 'bo')
axis equal
axis([-4 4 -4 4])
subplot(1,2,2)
plot(real(sample2), imag(sample2), 'ro')
axis equal
axis([-4 4 -4 4])
%% Testing the PDF
clc,clear
figure(3)
x = linspace(-5,5);
y1 = ComplexGaussianPDF(x, 0, 1);
y2 = CompoundGaussianPDF(x, 1);
hold on
plot(x, y1, 'b', LineWidth=1.5)
plot(x, y2, 'r', LineWidth=1.5)
legend('CN', 'CG')

% Guessing there's a normalisation factor missing in compound


%%
clc

sigma = 1;
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











