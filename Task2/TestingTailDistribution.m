clc
clear

x = linspace(0.1,10,1000);


sigma = 1;
N = 2;
nu = 1;


y1 = H_nKdist(x.^2, 0, sigma, nu);
y2 = H_nKdist(x.^2, 1, sigma, nu);
y3 = H_nKdist(x.^2, 5, sigma, nu);
y4 = H_nKdist(x.^2, 128, sigma, nu);
yG = H_nGaussian(x.^2, 0, sigma);

hold on
plot(x,y1, LineWidth = 1.5)
plot(x,y2, LineWidth = 1.5)
plot(x,y3, LineWidth = 1.5)
plot(x,y4, LineWidth = 1.5)
plot(x, yG, 'k--')
axis([0,10,0,10])
legend('K: h_0(x^2)', 'K: h_1(x^2)', 'h_5(x^2)','K: h_{128}(x^2)','CN: h_0(x^2)',location = 'best')
xlabel('x')





