function [U,cpu_time] = strangv(D,cl,gamma_par,kappa,zeta,U0,m,tau,m_k,xi)
%
% function [U,cpu_time] = strangv(D,cl,gamma_par,kappa,zeta,U0,m,tau,m_k,xi)
%
% Implementation of STRANG-V
  t = 0;

  tic
  n = [length(D{1}), length(D{2})];
  pn = prod(n);
  for mu = 1:2
    a{mu} = D{mu}(:,1);
    v{mu} = -xi*a{mu};
    v{mu}(1) = v{mu}(1)+1/2;
    cc{mu} = [v{mu};0;v{mu}(n(mu):-1:2)];
  end

  t = zeros(2*n);
  t(:,1) = cc{1};
  t(1,:) = t(1,:)+cc{2}';
  tev = fft2(t);

  % For tau preconditioner
  tau1{1} = D{1}(:,1)-[D{1}(3:n(1),1);0;0];
  tau1{2} = D{2}(:,1)-[D{2}(3:n(2),1);0;0];
  tauev = gentauev2(tau1);
  tauev1 = 1-xi*tauev;

  Mact = @(x) tx(tev,x);
  Mpre = @(x) tauinvx2(tauev1,x); % tau preconditioner

  U = U0(:);
  itvec = NaN(1,m);
  counter = 1;
  for jj = 1:m
    U = exp(gamma_par*(tau/2)-((kappa+1i*zeta)/(2*kappa))*log(1+2*(tau/2)*phiscal(2*gamma_par*(tau/2),1)*kappa*(real(U).^2+imag(U).^2))).*U;
    [U,mit] = phisil(tau*cl,Mact,Mpre,n,U,0,m_k,xi);
    itvec(counter) = mit;
    counter = counter + 1;
    U = exp(gamma_par*(tau/2)-((kappa+1i*zeta)/(2*kappa))*log(1+2*(tau/2)*phiscal(2*gamma_par*(tau/2),1)*kappa*(real(U).^2+imag(U).^2))).*U;
  end
  cpu_time = toc;
  disp(sprintf('Average iterations PCG: %.4f',mean(itvec)))

