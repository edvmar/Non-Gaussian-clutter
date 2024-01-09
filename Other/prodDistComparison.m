clear
clc


zValues = linspace(-3,3,500);

% complex Gaussian
sigma = 0.7;
mu    = 0;
complexGaussian = @(z) 1/(pi*sigma^2).*exp(-abs(z-mu).^2/sigma^2);

% complex Weibull
a  = 0.15; %sqrt(2);
nu_W = 0.8; %1/2;
complexWeibull = @(z) a*nu_W*abs(z).^(2*nu_W-2).*exp(-a*abs(z).^(2*nu_W));

% complex K
nu_K  = 1;
gamma = 1; %0.01;
complexK = @(z) 2*(nu_K/gamma)/(gamma(nu_K)*pi)*sqrt((nu_K/gamma)*abs(z).^2).^(nu_K-1)...
                .*besselk(nu_K-1,2*sqrt(nu_K/gamma*abs(z)));

plot(zValues,complexGaussian(zValues),zValues,complexWeibull(zValues), zValues, complexK(zValues))
xl = xlabel('x');
fontsize(xl,16,'points')
lgd = legend('complex Gaussian distribution', 'complex Weibull distribution', 'complex K distribution');
fontsize(lgd,16,'points')




