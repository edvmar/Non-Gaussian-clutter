%%

clc;clear;clf

 

N = 128;

ii = 0:(N-1);

d = 0.5/N^2;

v = exp(-d*ii.^2);

%S = toeplitz(v)+1e-10*eye(N);

S = CalculateToeplitzMatrix(N, d)+1e-10*eye(N);

imagesc(S)

colorbar()



%%
L = chol(S, 'lower');

Li = inv(L);

 

fprintf('\n\t%0.1e,\t %0.1e\n',det(S), det(L))

 

Si = inv(S);

SLi = Li'*Li;

fprintf('\n\t%0.1e,\t %0.1e,\t %0.1e\n',det(Si), det(SLi), det(Li))

 

%%

x = (randn(N,1)+1i*randn(N,1))/sqrt(2);

y = L*x;

 

s = 10*exp(1i*0.13*ii/N).';

 

Lis = (L\s);

q0  = x'*x

q1a = q0-2*real(x'*Lis)+Lis'*Lis

q1b = (L\(y-s)); q1b = real(q1b'*q1b)

q1c = real((y-s)'*SLi*(y-s))

q1d = real((y-s)'*Si*(y-s))