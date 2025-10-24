function T = tucker3d(T,A)
%
% function T = tucker3d(T,A)
%
% Implementation of 3D Tucker operator with square matrices using
% loop-over-GEMMs approach
  sT = size(T);
  T = reshape(A{1}*reshape(T,sT(1),[]),sT);
%  T = twomp3d(T,'none',A{2},'transpose'); % GNU Octave <= 10.2.0 and
                                          % MATLAB <= R2020a
  T = pagemtimes(T,'none',A{2},'transpose'); % MATLAB >= R2020b
  T = reshape(reshape(T,[],sT(3))*A{3}.',sT);

end
