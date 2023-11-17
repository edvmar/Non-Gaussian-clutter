clear
clc

%% Testing the sample

sample = SampleCompoundGaussian(10^3, 1, 1);

figure(1)
histogram(abs(sample))

%% Testing the PDF

figure(2)
x = linspace(0,10);
y = CompoundGaussianPDF(x, 1, 1);
plot(x,y)