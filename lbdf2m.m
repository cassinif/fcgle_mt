function [U,cpu_time] = lbdf2m(D,cl,g,U0,m,tau)
%
% function [U,cpu_time] = lbdf2m(D,cl,g,U0,m,tau)
%
% Implementation of LBDF2-M
  t = 0;

  tic
  [Q{1},L{1}] = eig(D{1},'vector');
  [Q{2},L{2}] = eig(D{2},'vector');
  [LL{1},LL{2}]=ndgrid(L{1},L{2});
  Mm = LL{1}+LL{2};
  Mm1 = (1-((2*tau/3)*(cl*Mm)));

  Up = U0;

  Utmp = Q{1}*(((Q{1}.'*((U0+tau*g(0,U0))*Q{2}))./(1-(tau*(cl*Mm))))*Q{2}.');
  t = t + tau;

  for jj = 2:m
    gn = g(t+tau,2*Utmp-Up);

    U = Q{1}*(((Q{1}.'*(((4/3)*Utmp-Up/3+(2*tau/3)*gn)*Q{2}))./Mm1)*Q{2}.');

    Up = Utmp;
    Utmp = U;
    t = t + tau;
  end
  cpu_time = toc;
