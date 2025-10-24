function Y = dstI2(B)
%
% function Y = dstI2(B)
%
% two-dimensional DST-I of matrix B
n = size(B);
% Odd extension of B
Bext = [zeros(1,2*(n(2)+1));...
        zeros(n(1),1),B,zeros(n(1),1),-B(:,n(2):-1:1);...
        zeros(1,2*(n(2)+1));...
        zeros(n(1),1),-B(n(1):-1:1,:),zeros(n(1),1),B(n(1):-1:1,n(2):-1:1)];
% FFT
Y = fft2(Bext);
% Extract DST-I values (real parts)
Y = -Y(2:n(1)+1,2:n(2)+1)/4;
