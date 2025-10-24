function ret = phiscal(z,ell)
%
% function ret = phiscal(z,ell)
%
% This function returns the scalar phi_\ell function computed either by Taylor
% series or by the recurrence formula
  if ell == 0
    ret = exp(z);
  else
    idx = abs(z) < 1;
    ret = zeros(size(z));
    ret(idx) = polyval(1./factorial(16+ell:-1:ell),z(idx));
    ret(~idx) = exp(z(~idx)) - 1;
    fact = 1;
    for jj = 1:(ell-1)
      fact = fact*jj;
      ret(~idx) = ret(~idx) - z(~idx).^jj/fact;
    end
    ret(~idx) = ret(~idx)./(z(~idx).^ell);
  end
