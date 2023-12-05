
nus = [1, 2, 3, 50];

x = linspace(0,3);
figure(1)
hold on

for i = 1:length(nus)

    nu = nus(i);
    y = besselk(nu,x);
    y = TailDistributionCompoundGaussian(x, 0, 1, nu);

    plot(x,y, LineWidth=1)
end

y2 = TailDistributionComplexGaussian(x, 0, 1);
plot(x,y2,'k--',LineWidth=1)

axis([0, x(100), 0, 1])
legend('nu = 1','nu = 2','nu = 3','nu = 50', 'Gauss')


% Seems to not be defined for nu > 171. Approaches Gaussian as nu ≥ 30.

