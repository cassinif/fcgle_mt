function ev = genl2ev(t)
%
% function ev = genl2ev(t)
%
% genl2ev computes eigenvalues of the BCCB preconditioner.
[n,m] = size(t);
n = n/2;
v = zeros(n,m);
v(1,:) = t(1,:);
i = (2:n)';
v(i,:) = ((n-(i-1)).*t(i,:)+(i-1).*t(n+i,:))/n;
v = fft(v); % CB eigenvalues.
v = v.';
m = m/2;
ev(1,:) = v(1,:);
i = (2:m)';
ev(i,:) = ((m-(i-1)).*v(i,:)+(i-1).*v(m+i,:))/m;
if min(n,m) >= 2
  ev = fft(ev);
end
ev = ev.';
