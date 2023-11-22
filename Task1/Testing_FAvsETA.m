clear
clc

sampleSize = 1e8;
sample = SampleComplexGaussian(sampleSize, 0, 1);

SIR = 13;
SIR = 10^(SIR/10);           % potentially like this ? 
alpha = sqrt(SIR);             
   
theta = 0; 
s = alpha*(cos(theta)+1i*sin(theta)); % signal 

% [d,index] = min(abs(sample-s));
% [D,Index] = max(abs(sample-s));
% 
% x = sample(index)
% X = sample(Index)
% 
% fH0 = ComplexGaussianPDF(x,0,1);
% fH1 = ComplexGaussianPDF(x,s,1);
% LRTmin = fH1/fH0
% 
% FH0 = ComplexGaussianPDF(X,0,1);
% FH1 = ComplexGaussianPDF(X,s,1);
% LRTmax = FH1/FH0


eta = 4*1e5;

a_l = (log(eta)+alpha^2)/(2*alpha);


fH1_fa = ComplexGaussianPDF(sample, s, 1);           
fH0_fa = ComplexGaussianPDF(sample, 0, 1);
sumFAnaiv = 0;
sumFAthresh = 0;
a_l_numerical = 100;
for i = 1:sampleSize
    if fH1_fa(i)/fH0_fa(i) > eta
        %disp(['x = ', num2str(sample(i)), ' fH1 = ', num2str(fH1_fa(i)), ' fH0 = ', num2str(fH0_fa(i)), ' LRT = ', num2str(fH1_fa(i)/fH0_fa(i))])
        sumFAnaiv = sumFAnaiv + 1; 
        if (a_l_numerical > real(sample(i)))
            a_l_numerical = real(sample(i));
        end
    end
    if (real(sample(i)) > a_l)
        sumFAthresh = sumFAthresh + 1;
    end
end
sumFA = sum(((fH1_fa./fH0_fa) > eta));
sumFA
sumFAnaiv
sumFAthresh

Numerical = sumFA/sampleSize
Numerical_WithThresh = sumFAthresh/sampleSize
Analytical = 1 - normcdf(sqrt(2)*a_l)



%%
hold on
plot(real(sample+s), imag(sample+s), 'ro')
plot(real(sample), imag(sample), 'bo')
plot([a_l, a_l], [-5,5],'k--', LineWidth = 1.0)
plot([a_l_numerical, a_l_numerical], [-5,5],'m--', LineWidth = 1.0)

