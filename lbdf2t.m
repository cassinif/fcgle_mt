function [U,cpu_time] = lbdf2t(D,cl,g,U0,m,tau)
%
% function [U,cpu_time] = lbdf2t(D,cl,g,U0,m,tau)
%
% Implementation of LBDF2-T
  t = 0;

  tic
  [Q{1},L{1}] = eig(D{1},'vector');
  [Q{2},L{2}] = eig(D{2},'vector');
  [Q{3},L{3}] = eig(D{3},'vector');
  [LL{1},LL{2},LL{3}]=ndgrid(L{1},L{2},L{3});
  Mm = LL{1}+LL{2}+LL{3};
  Mm1 = (1-((2*tau/3)*(cl*Mm)));

  Up = U0;

  Utmp = tucker3d(ttucker3d(U0 + tau*g(0,U0),Q)./(1-(tau*(cl*Mm))),Q);
  t = t + tau;

  for jj = 2:m
    gn = g(t+tau,2*Utmp-Up);

    U = tucker3d(ttucker3d((4/3)*Utmp-(1/3)*Up+(2*tau/3)*gn,Q)./Mm1,Q);

    Up = Utmp;
    Utmp = U;
    t = t + tau;
  end
  cpu_time = toc;
