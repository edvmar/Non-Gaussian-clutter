

clc
clear

x = linspace(0.01,5);


sigma = 1;
N = 12;
nu = 1;


y = H_nKdist(x.^2, N, sigma, nu);
y2 = H_nGaussian(x.^2, 0, sigma);

hold on
plot(x,y, LineWidth = 1.5)
plot(x, y2, 'k--')
axis([0,5,0,2])




