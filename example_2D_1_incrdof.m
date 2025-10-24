clear all
close all

% Example 2D number 1: TDSRC with LBDF2 using the proposed technique
%                      or GMRES + tau precond
% Increasing space DOF, fixed number of steps

d = 2;

nrange = 200:100:600;

alpha = [1.2,1.8];

nsteps = 25;

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

tau = T/nsteps;

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

  S = @(t) tmpten*exp(-1i*t);
  G = @(t,u) gamma_par*u -coeff_nonlin.*(real(u).^2+imag(u).^2).*u + S(t);

  s = @(t) tmpten(:)*exp(-1i*t);
  g = @(t,u) gamma_par*u -coeff_nonlin.*(real(u).^2+imag(u).^2).*u + s(t);

  disp('LBDF2 M')
  [Ulbdf2m, cpu_lbdf2m(counter)] = lbdf2m(D,coeff_lin,G,U0,nsteps,tau);

  disp('LBDF2 V')
  [Ulbdf2v,cpu_lbdf2v(counter)] = lbdf2v(D,coeff_lin,g,U0(:),nsteps,tau);
end

figure;
loglog(nrange,cpu_lbdf2v,'-or')
hold on
loglog(nrange,cpu_lbdf2m,'-+b')
legend('LBDF2 V','LBDF2 M')
xlabel('n')
ylabel('Wall-clock time')
title('Increasing number of dof')
drawnow
