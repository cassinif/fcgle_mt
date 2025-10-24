clear all
close all

% Example 3D number 2: NOSRC with STRANG using the proposed technique
%                      in SINGLE precision
% Fixed space DOF, increasing number of steps

d = 3;

n = 250*ones(1,d);

alpha = [1.2,1.8,1.5];

mrange = 5:5:25;

a = -10*ones(1,d);
b = -a;

nu = 1;
eta = 1;
gamma_par = 1;
kappa = 1;
zeta = 1;
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
  D{mu} = toeplitz(gg);
end
[X{1:d}] = ndgrid(x{:});

U0 = sech(X{1}).*exp(1i*X{1});
for mu = 2:d
  U0 = U0.*sech(X{mu}).*exp(1i*X{mu});
end

g = @(t,u) gamma_par*u -coeff_nonlin.*(real(u).^2+imag(u).^2).*u;


nref = 300;
disp('Computing reference solution...')
Uref = strangt(D,coeff_lin,gamma_par,kappa,zeta,g,U0,nref,T/nref);
disp('Reference solution computed!')

coeff_nonlin = single(coeff_nonlin);
gamma_par = single(gamma_par);
g = @(t,u) gamma_par*u -coeff_nonlin.*(real(u).^2+imag(u).^2).*u;

for mu = 1:d
  D{mu} = single(D{mu});
end
coeff_lin = single(coeff_lin);
kappa = single(kappa);
zeta = single(zeta);
U0 = single(U0);

counter = 0;
for nsteps = mrange
  nsteps
  counter = counter + 1;
  tau = T/nsteps;

  disp('STRANG T')
  [Ust, cpu_st(counter)] = strangt(D,coeff_lin,gamma_par,kappa,zeta,g,U0,nsteps,single(tau));
  err_st(counter) = sqrt(prod(h)*sum(abs(Ust(:)-Uref(:)).^2));

end

figure;
loglog(err_st,cpu_st,'-+b')
legend('STRANG T')
xlabel('Error')
ylabel('Wall-clock time')
title('Work-precision diagram')
drawnow
