function x = idstI2(y)
% Inverse 2D DST-I (matches the unnormalized forward)
% Uses the fact that S is symmetric and S*S' = (M+1)/2 * I in 1D.
% So X = (4/((M+1)(N+1))) * DST-I(Y).

    n = size(y);
    x = 4/prod(n+1)*dstI2(y);
end
