%% Uppgift 1a)
%% Gaussian and Gaussian
clc,clear

n = 10^4; % Number of samples
alpha = 1;
theta = 0;
s = alpha*(cos(theta)+1i*sin(theta));
sigma_detector = 1; % The standard deviation for the detector
sigma_clutter = 1; % The standard deviation for the detector
mean_detector = 0;
mean_clutter = 0;

samples = randn(1,n); % ?
a = samples*sigma_clutter + mean_clutter;
samples = randn(1,n);
b = samples*sigma_clutter + mean_clutter;
c = a + 1i*b;
c_s = c + s;

number_of_eta_values = 20;
eta_values = linspace(1,10,number_of_eta_values);
Pfa = zeros(1, number_of_eta_values);
Pda = zeros(1, number_of_eta_values);

G_pdf = @(z, mu, sigma) 1/(pi*sigma^2)*exp(-abs(z-mu)^2/sigma^2); % Is this correct?


for i=1:length(eta_values)
    cur_eta = eta_values(i);
    sum_fa = 0;
    sum_pda = 0;
    for j=1:n
        if(G_pdf(c(j), s, sigma_detector)/G_pdf(c(j), 0, sigma_detector)>cur_eta)
            sum_fa = sum_fa + 1;
        end

        if(G_pdf(c_s(j), s, sigma_detector)/G_pdf(c_s(j), 0, sigma_detector)>cur_eta)
            sum_pda = sum_pda + 1;
        end
    end
    Pfa(i) = sum_fa/n;
    Pda(i) = sum_pda/n;
end
plot(Pfa,Pda)
%% Gaussian and Compound Gaussian

clc,clear

n = 10^4; % Number of samples
alpha = 1;
theta = 0;
s = alpha*(cos(theta)+1i*sin(theta));
sigma_detector = 1; % The standard deviation for the detector
sigma_clutter = 1; % The standard deviation for the detector
mean_detector = 0;
mean_clutter = 0;


pdf = @(x,v) v/(gamma(v)*2^(v-2)*pi)*(2*sqrt(v)*abs(x)).^(v-1).*bessely(v-1,2*sqrt(v)*abs(x));
x=linspace(-1,1,100);
plot(x,pdf(x, 0.50))
%%

pdf = @(x, b, v) 2*b^v/gamma(v) * x.^(v-1) .* exp(-b * x) .* bessely(v-1, 2*sqrt(b*x));

% Example usage:
b = 2;
v = 3;
x_values = linspace(0, 10, 1000);

pdf_values = pdf(x_values, b, v);

figure;
plot(x_values, pdf_values, 'LineWidth', 2);
title('K Distribution PDF');
xlabel('X');
ylabel('Probability Density');
grid on;

%%

samples = randn(1,n); % ?
a = samples*sigma_clutter + mean_clutter;
samples = randn(1,n);
b = samples*sigma_clutter + mean_clutter;
c = a + 1i*b;
c_s = c + s;

number_of_eta_values = 20;
eta_values = linspace(1,10,number_of_eta_values);
Pfa = zeros(1, number_of_eta_values);
Pda = zeros(1, number_of_eta_values);

G_pdf = @(z, mu, sigma) 1/(pi*sigma^2)*exp(-abs(z-mu)^2/sigma^2); % Is this correct?


for i=1:length(eta_values)
    cur_eta = eta_values(i);
    sum_fa = 0;
    sum_pda = 0;
    for j=1:n
        if(G_pdf(c(j), s, sigma_detector)/G_pdf(c(j), 0, sigma_detector)>cur_eta)
            sum_fa = sum_fa + 1;
        end

        if(G_pdf(c_s(j), s, sigma_detector)/G_pdf(c_s(j), 0, sigma_detector)>cur_eta)
            sum_pda = sum_pda + 1;
        end
    end
    Pfa(i) = sum_fa/n;
    Pda(i) = sum_pda/n;
end
plot(Pfa,Pda)

%% Uppgift 1a analytiskt
a = (eta+alpha^2)/(2*alpha)
% Vi ska integrera från a till inf (real), och -inf till inf (imaginär del). 

\int e^(-x^2-y^2)/pi dxdy
\int sqrt(2pi)/pi e^(-(x)^2)
\int sqrt(2pi)/pi e^(-(x-alpha)^2)
2*K_o(2sqrt(x))



F^-1(0.5) = y
0.5 = F(y)
int(-inf,y) (k_distribution) = 0.5





