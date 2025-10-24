function [u,cpu_time] = lbdf2v(D,cl,g,u0,m,tau)
%
% function [u,cpu_time] = lbdf2v(D,cl,g,u0,m,tau)
%
% Implementation of LBDF2-V
  t = 0;

  tic
  n(1) = length(D{1});
  n(2) = length(D{2});
  t1 = -cl*[D{1}(:,1);0;D{1}(1,n(1):-1:2)'];
  t2 = -cl*[D{2}(1,:),0,D{2}(n(2):-1:2,1)'];
  tt = zeros(2*n(1),2*n(2));
  tt(:,1) = tau*t1;
  tt(1,:) = tt(1,:)+tau*t2;
  tt(1,1) = 1+tt(1,1);
  tev = fft2(tt);
  Tx = @(x) tx(tev,x); % BTTB action

  g0 = g(t,u0);

  tau1{1} = D{1}(:,1)-[D{1}(3:n(1),1);0;0];
  tau1{2} = D{2}(:,1)-[D{2}(3:n(2),1);0;0];
  tauev = gentauev2(tau1);
  tauev1 = 1-(tau*cl)*tauev;
  tauinvx = @(x) tauinvx2(tauev1,x); % tau preconditioner

  itvec = NaN(1,m);
  ct = 1;

  % First step
  [utmp,flag,relres,iter] = gmres(Tx,u0(:)+tau*g0,20,[],1,tauinvx,[],u0);
  itvec(ct) = iter(2);
  ct = ct + 1;

  t = t + tau;

  tt = zeros(2*n(1),2*n(2));
  tt(:,1) = (2*tau/3)*t1;
  tt(1,:) = tt(1,:)+(2*tau/3)*t2;
  tt(1,1) = 1+tt(1,1);
  tev = fft2(tt);

  Tx = @(x) tx(tev,x); % BTTB action

  tauevn = 1-(2*tau/3*cl)*tauev;
  tauinvx = @(x) tauinvx2(tauevn,x); % tau preconditioner

  up = u0;
  for jj = 2:m
    gn = g(t+tau,2*utmp-up);
    [u,flag,relres,iter] = gmres(Tx,4/3*utmp-up/3+(2*tau/3)*gn,20,[],1,tauinvx,[],utmp);
    itvec(ct) = iter(2);
    ct = ct + 1;
    up = utmp;
    utmp = u;
    t = t + tau;
  end
  cpu_time = toc;
  disp(sprintf('Average iterations GMRES: %.4f',mean(itvec)))
