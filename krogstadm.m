function [U,cpu_time] = krogstadm(D,cl,g,U0,m,tau)
%
% function [U,cpu_time] = krogstadm(D,cl,g,U0,m,tau)
%
% Implementation of KROGSTAD-M
  t = 0;

  tic
  [Q{1},L{1}] = eig(D{1},'vector');
  [Q{2},L{2}] = eig(D{2},'vector');
  [LL{1},LL{2}]=ndgrid(L{1},L{2});
  Mm = LL{1}+LL{2};

  P3 = phiscal((tau*cl)*Mm,3);
  P2 = phiscal((tau*cl)*Mm,2);
  P2h = phiscal((tau/2*cl)*Mm,2);
  P1 = phiscal((tau*cl)*Mm,1);
  P1h = phiscal((tau/2*cl)*Mm,1);

  U = U0;
  for jj = 1:m
    gn = g(t,U);
    ttfn = Q{1}.'*((cl*(D{1}*U+U*D{2}.')+gn)*Q{2});
    U2 = U + tau/2*(Q{1}*((P1h.*ttfn)*Q{2}.'));
    ttgn2gn = Q{1}.'*((g(t+tau/2,U2)-gn)*Q{2});
    U3 = U2 + tau*(Q{1}*((P2h.*ttgn2gn)*Q{2}.'));
    ttgn3gn = Q{1}.'*((g(t+tau/2,U3)-gn)*Q{2});
    U4 = U + tau*(Q{1}*((P1.*ttfn+(2*P2).*ttgn3gn)*Q{2}.'));
    ttgn4gn = Q{1}.'*((g(t+tau,U4)-gn)*Q{2});
    U = U + tau*(Q{1}*((P1.*ttfn+(2*P2-4*P3).*ttgn2gn+...
        (2*P2-4*P3).*ttgn3gn+(-P2+4*P3).*ttgn4gn)*Q{2}.'));
    t = t + tau;
  end
 cpu_time = toc;
