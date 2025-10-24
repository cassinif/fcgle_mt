function y = tauinvx2(tauev,d)
%
% function y = tauinvx2(tauev,d)
%
% solves the preconditioning (tau) system tau y = d, tau the Kronecker
% sum of two preconditioners
% tauev: the eigenvalues computed by gentauev2
[n,m] = size(tauev);
y = reshape(idstI2(dstI2(reshape(d,n,m))./tauev),n*m,1);
