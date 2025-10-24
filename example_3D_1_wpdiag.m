clear all
close all

% Example 3D number 1: TDSRC with LBDF2 using the proposed technique
%                      in SINGLE precision
% Fixed space DOF, increasing number of steps

d = 3;

n = 200*ones(1,d);

alpha = [1.2,1.8,1.5];

mrange = 15:5:35;

a = -1*ones(1,d);
b = -a;

nu = 1;
eta = 1;
gamma_par = 3;
kappa = 1;
zeta = 2;
T = 1;

coeff_lin = nu+1i*eta;
coeff_nonlin = kappa+1i*zeta;

for mu = 1:d
  x{mu} = linspace(a(mu),b(mu),n(mu)+2).';
  x{mu} = x{mu}(2:n(mu)+1);
  h(mu) = (b(mu)-a(mu))/(n(mu)+1);
  ii = (1:n(mu)-2)';
  gg = -1/h(mu)^alpha(mu)*gamma(alpha(mu)+1)./...
        [gamma(alpha(mu)/2+1)^2;-gamma(alpha(mu)/2)*gamma(alpha(mu)/2+2);...
        (-1).^(2:n(mu)-1)'*gamma(alpha(mu)/2)*gamma(alpha(mu)/2+2).*...
        cumprod((alpha(mu)/2+1+ii)./(alpha(mu)/2-ii))];
  D{mu} = single(toeplitz(gg));
end
[X{1:d}] = ndgrid(x{:});

U0 = ones(n);
for mu = 1:d
  U0 = U0.*(X{mu}.^2+a(mu)).^4;
end

tmpten = (-1i-gamma_par+coeff_nonlin*(U0.^2)).*U0;

for mu = 1:d
  tmp = coeff_lin*ones(n);
  for kk = [1:mu-1,mu+1:d]
    tmp = tmp.*(X{kk}.^2+a(kk)).^4;
  end
  tmp2 = zeros(n);
  cc = [1,-8,24,-32,16];
  ct = 0;
  for ll = 9:-1:5
    ct = ct + 1;
    tmp2 = tmp2 + cc(ct)*gamma(ll)/gamma(ll-alpha(mu))*((b(mu)+X{mu}).^(ll-1-alpha(mu))+(b(mu)-X{mu}).^(ll-1-alpha(mu)));
  end
  tmpten = tmpten + (tmp/(2*cos(alpha(mu)*pi/2))).*tmp2;
end

tmpten = single(tmpten);
s = @(t) tmpten*exp(-1i*t);
coeff_nonlin = single(coeff_nonlin);
gamma_par = single(gamma_par);
g = @(t,u) gamma_par*u -coeff_nonlin.*(real(u).^2+imag(u).^2).*u + s(t);

Uref = U0.*exp(-1i*T);

coeff_lin = single(coeff_lin);
U0 = single(U0);

counter = 0;
for nsteps = mrange
  nsteps
  counter = counter + 1;
  tau = T/nsteps;

  disp('LBDF2 T')
  [Ulbdf2t, cpu_lbdf2t(counter)] = lbdf2t(D,coeff_lin,g,U0,nsteps,single(tau));
  err_lbdf2t(counter) = sqrt(prod(h)*sum(abs(Ulbdf2t(:)-Uref(:)).^2));
end

figure;
loglog(err_lbdf2t,cpu_lbdf2t,'-xb')
legend('LBDF2 T')
xlabel('Error')
ylabel('Wall-clock time')
title('Work-precision diagram')
drawnow
