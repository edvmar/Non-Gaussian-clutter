%%%%%%%%%%%%% Sample Testing %%%%%%%%%%%%%%%%%
%
% Trying to implement the sampling 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = @(x) 1 - exp(-abs(x).^2);  % eqn (12)

radii = linspace(0, 5, 1000);
values = F(radii);

n = 1e4;
uniforms = rand(1,n);
sample = [];

for i = 1:n
    uniform = uniforms(i);
    [~, index] = min(abs(values-uniform));
    radius = radii(index);
    
    noPoints = fix(radius);
    rest = radius - fix(radius);
    if (rand(1,1) < rest)
        noPoints = noPoints + 1;
    end
    %noPoints = 1;
    thetas = rand(1,noPoints)*2*pi;
    sample = [sample, exp(thetas*1i)*radius];
    
end
figure(1)
histogram2(real(sample), imag(sample))
figure(2)
plot(real(sample), imag(sample), 'ro')
%%
figure(3)
hold on
histogram(real(sample), 'Normalization','pdf')
x = linspace(-4,4);
sigma = 1;
f = @(x) 1/(sqrt(2*pi*sigma^2))*exp(-x.^2/sigma^2);
plot(x,f(x))

