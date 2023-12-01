%%%%%%%%%%%%% Sample Testing %%%%%%%%%%%%%%%%%
%
% Trying to implement the sampling 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 1;
F = @(x) 1 - exp(-abs(x).^2/sigma^2);  % eqn (12)

xmax = 4*sigma;
radii = linspace(0, xmax, 1000);
values = F(radii);

n = 1e3;
uniforms = rand(1,n);
sample1 = zeros(1,n);

for i = 1:n
    uniform = uniforms(i);
    [~, index] = min(abs(values-uniform));
    radius = radii(index);
    
    noPoints = 1;
    thetas = rand(1,noPoints)*2*pi;
    sample1(i) = exp(thetas*1i)*radius; 
end

%% K - dist
nu = 1;
eta = sigma^2;

F = @(x) 1 - (2*(sqrt(nu/eta).*abs(x)).^nu)/gamma(nu).*besselk(nu,2*sqrt(nu/eta)*abs(x));  % eqn (12)   (maybe nu-1 or nu in Bessel ??)

radii = linspace(0, 2*xmax, 1000);
values = F(radii);

uniforms = rand(1,n);
sample2 = zeros(1,n);

for i = 1:n
    uniform = uniforms(i);
    [~, index] = min(abs(values-uniform));
    radius = radii(index);
    
    noPoints = 1;
    thetas = rand(1,noPoints)*2*pi;
    sample2(i) = exp(thetas*1i)*radius; 
end

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
%xmax=2*xmax;
axis([-xmax xmax -xmax xmax])
subplot(1,2,2)
plot(real(sample2), imag(sample2), 'ro')
axis equal
axis([-xmax xmax -xmax xmax])
%
x = linspace(-xmax,xmax);
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
