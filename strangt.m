function [U,cpu_time] = strangt(D,cl,gamma_par,kappa,zeta,gnl,U0,m,tau)
%
% function [U,cpu_time] = strangt(D,cl,gamma_par,kappa,zeta,gnl,U0,m,tau)
%
% Implementation of STRANG-T
  t = 0;

  tic
  [Q,L] = eig(D{1},'vector');
  E{1} = Q*(exp((tau*cl)*L).*Q.');
  [Q,L] = eig(D{2},'vector');
  E{2} = Q*(exp((tau*cl)*L).*Q.');
  [Q,L] = eig(D{3},'vector');
  E{3} = Q*(exp((tau*cl)*L).*Q.');

  U = U0;
  for jj = 1:m
    U = exp(gamma_par*(tau/2)-((kappa+1i*zeta)/(2*kappa))*log(1+2*(tau/2)*phiscal(2*gamma_par*(tau/2),1)*kappa*(real(U).^2+imag(U).^2))).*U;
    U = tucker3d(U,E);
    U = exp(gamma_par*(tau/2)-((kappa+1i*zeta)/(2*kappa))*log(1+2*(tau/2)*phiscal(2*gamma_par*(tau/2),1)*kappa*(real(U).^2+imag(U).^2))).*U;
  end
  cpu_time = toc;
