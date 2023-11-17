%Matthew Ritchie
%16/06/10

%K-distribution pfa
function [pfadB,nu] = RDV_kdist6(RDV_Pfa_X_pos)
%v = shape parameter
%Et = Threshold
%b = Scale parameter
%Kv = Modified Bessel Function
%R = Gamma function
a = RDV_Pfa_X_pos;
nu = [0.1:0.1:10,100];
pfadB = zeros(length(a),length(nu));
%selecting value of nu from array
v = nu;
G = gamma(v);
a = 10.^(a./10);
[aa,GG]=meshgrid(a,G);
%Modified Bessel function of 2nd kind
[aa,vv] = meshgrid(a,v);
H = 2*((aa.*vv).^0.5);
[Kv] = besselk(vv,H);
pfa = (2./GG).*((aa.*vv).^(vv./2)).*Kv;
pfadB(:,:) = log10(pfa)';





