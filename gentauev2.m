function tauev = gentauev2(t)
%
% function gentauev2(t)
%
% gentauev2 computes the eigenvalues of the Kronecker sum of
% two tau preconditioners
d = 2;
for mu = 1:d
  n(mu) = length(t{mu});
  alpha{mu} = pi*(1:n(mu))'/(n(mu)+1);
  tauev{mu} = sum(t{mu}'.*sin(alpha{mu}*(1:n(mu))),2)./sin(alpha{mu});
end
tauev = tauev{1}+tauev{2}.';
