function [U,cpu_time] = krogstadv(D,cl,g,U0,m,tau,m_k,xi)
%
% function [U,cpu_time] = krogstadv(D,cl,g,U0,m,tau,m_k,xi)
%
% Implementation of KROGSTAD-V
  t = 0;

  tic
  n = [length(D{1}), length(D{2})];
  pn = prod(n);

  for mu = 1:2
    a{mu} = D{mu}(:,1);
    cca{mu} = [a{mu};0;a{mu}(n(mu):-1:2)];
    v{mu} = -xi*a{mu};
    v{mu}(1) = v{mu}(1)+1/2;
    cc{mu} = [v{mu};0;v{mu}(n(mu):-1:2)];
  end
  tt = zeros(2*n);
  tt(:,1) = cca{1};
  tt(1,:) = tt(1,:)+cca{2}';
  teva = fft2(tt);

  tt = zeros(2*n);
  tt(:,1) = cc{1};
  tt(1,:) = tt(1,:)+cc{2}';
  tev = fft2(tt);
  ev2 = genl2ev(tt);

  % For tau preconditioner
  tau1{1} = D{1}(:,1)-[D{1}(3:n(1),1);0;0];
  tau1{2} = D{2}(:,1)-[D{2}(3:n(2),1);0;0];
  tauev = gentauev2(tau1);
  tauev1 = 1-xi*tauev;

  U = U0;

  Mact = @(x) tx(tev,x);
  Mpre = @(x) tauinvx2(tauev1,x); % tau preconditioner

  mit = NaN(5,m);
  ct = 1;
  for jj = 1:m
    gn = g(t,U);
    fn = cl*tx(teva,U) + gn;
    [phifn,mit(1,ct)] = phisil([tau/2,tau]*cl,Mact,Mpre,n,fn,1,m_k,xi);
    U2 = U + tau/2*phifn(:,1);
    D2 = g(t+tau/2,U2) - gn;
    [tmp,mit(2,ct)] = phisil(tau/2*cl,Mact,Mpre,n,D2,2,m_k,xi);
    U3 = U2 + tau*tmp;
    D3 = g(t+tau/2,U3) - gn;
    [tmp,mit(3,ct)] = phisil(tau*cl,Mact,Mpre,n,D3,2,m_k,xi);
    U4 = U + tau*phifn(:,2) + (2*tau)*tmp;
    D4 = g(t+tau,U4) - gn;

    [tmp1,mit(4,ct)] = phisil(tau*cl,Mact,Mpre,n,2*D2+2*D3-D4,2,m_k,xi);
    [tmp2,mit(5,ct)] = phisil(tau*cl,Mact,Mpre,n,-4*D2-4*D3+4*D4,3,m_k,xi);
    ct = ct + 1;
    U = U + tau*phifn(:,2) + tau*tmp1+tau*tmp2;
    t = t + tau;
  end
  cpu_time = toc;

  disp(sprintf('Average iterations PCG: %.4f',mean(mean(mit))))

end
