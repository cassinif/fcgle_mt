clear all
close all

% Example 2D number 2: NOSRC with STRANG using the proposed technique
%                      or SIL + PCG + tau precond
% Increasing DOF, fixed number of steps

d = 2;

nrange = 400:200:1200;

alpha = [1.2,1.8];

nsteps = 15;

m_k = 10;

a = -10*ones(1,d);
b = -a;

nu = 1;
eta = 1;
gamma_par = 1;
kappa = 1;
zeta = 1;
T = 1;

tau = T/nsteps;
xi = tau/10;

coeff_lin = nu+1i*eta;
coeff_nonlin = kappa+1i*zeta;

counter = 0;
for nn = nrange
  nn
  counter = counter + 1;
  n = nn*ones(1,d);
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

  disp('STRANG M')
  [Usm, cpu_sm(counter)] = strangm(D,coeff_lin,gamma_par,kappa,zeta,U0,nsteps,tau);

  disp('STRANG V')
  [Usv, cpu_sv(counter)] = strangv(D,coeff_lin,gamma_par,kappa,zeta,U0,nsteps,tau,m_k,xi);

end

figure;
loglog(nrange,cpu_sv,'-or')
hold on
loglog(nrange,cpu_sm,'-+b')
legend('STRANG V','STRANG M')
xlabel('n')
ylabel('Wall-clock time')
title('Increasing number of dof')
drawnow
