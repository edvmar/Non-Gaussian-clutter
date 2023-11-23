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
x = linspace(-10,10);
y1 = ComplexGaussianPDF(x, 0, 1);
y2 = CompoundGaussianPDF(x, 1, 1);
hold on
plot(x, y1, 'b', LineWidth=1.5)
plot(x, y2, 'r', LineWidth=1.5)
legend('CN', 'CG')

% Guessing there's a normalisation factor missing in compound

%%


