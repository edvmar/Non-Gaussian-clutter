%%%%%%%%%%%%% Sample Testing %%%%%%%%%%%%%%%%%
%
% Trying to implement the sampling 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = @(x) 1 - exp(-abs(x).^2);  % eqn (12)

radii = linspace(0, 5, 1000);
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
eta = 1;

F = @(x) 1 - (2*(sqrt(nu/eta).*abs(x)).^nu)/gamma(nu).*besselk(nu,2*sqrt(nu/eta)*abs(x));  % eqn (12)   (maybe nu-1 or nu in Bessel ??)

radii = linspace(0, 5, 1000);
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



figure(1)
subplot(1,2,1)
histogram2(real(sample1), imag(sample1))
subplot(1,2,2)
histogram2(real(sample2), imag(sample2))

figure(2)
subplot(1,2,1)
plot(real(sample1), imag(sample1), 'bo')
axis equal
axis([-4 4 -4 4])
subplot(1,2,2)
plot(real(sample2), imag(sample2), 'ro')
axis equal
axis([-4 4 -4 4])
%
x = linspace(-4,4);
sigma = 1/sqrt(2);
f = @(x) 1/(sqrt(2*pi*sigma^2))*exp(-x.^2/(2*sigma^2));
figure(3)

subplot(1,2,1)
hold on
histogram(real(sample1), 'Normalization','pdf')
plot(x,f(x))


subplot(1,2,2)
hold on
histogram(real(sample2), 'Normalization','pdf')
plot(x,f(x))
