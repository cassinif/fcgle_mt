function [out,mit] = phisil(taurange,Mact,Mpre,n,v,p,m,xi)
%
% function [out,mit] = phisil(taurange,Mact,Mpre,n,v,p,m,xi)
%
% This function returns the action of a phi_p function using the
% shift-and-invert Lanczos method. The inputs are coherent with the FCGLE
% under consideration

  ct = 1;
  itt = NaN(1,m);

  pn = prod(n);
  V = zeros(pn,m);
  H = zeros(m,m);

  z = zeros(pn,1); % initial guess for pcg

  beta = sqrt(v'*v);
  V(:,1) = v/beta;

  % j = 1
  [v,flag,relres,itt(ct)] = pcg(Mact,V(:,1),[],20,Mpre,[],z);
  ct = ct + 1;
  H(1,1) = V(:,1)'*v;
  v = v - H(1,1)*V(:,1);
  H(2,1) = sqrt(v'*v);
  H(1,2) = H(2,1);
  V(:,2) = v/H(2,1);

  for j = 2:(m-1)
    [v,flag,relres,itt(ct)] = pcg(Mact,V(:,j),[],20,Mpre,[],z);
    ct = ct + 1;

    v = v - H(j-1,j)*V(:,j-1);
    H(j,j) = V(:,j)'*v;
    v = v - H(j,j)*V(:,j);
    H(j+1,j) = sqrt(v'*v);
    H(j,j+1) = H(j+1,j);
    V(:,j+1) = v/H(j+1,j);
  end

  % j = m
  [v,flag,relres,itt(ct)] = pcg(Mact,V(:,m),[],20,Mpre,[],z);
  ct = ct + 1;
  v = v - H(m-1,m)*V(:,m-1);
  H(m,m) = V(:,m)'*v;

  [VV,DD] = eig(real(H),'vector'); % just spurious imaginary part in our case
  out = zeros(pn,length(taurange));
  for i = 1:length(taurange)
    out(:,i) = beta*(V*(VV*(phiscal((taurange(i)/xi)*(1-1./DD(:)),p).*reshape(VV(1,:),m,1))));
  end

  mit = mean(itt);
end
