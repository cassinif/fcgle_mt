function T=twomp3d(T,s,A,s2)
%
% function T=twomp3d(T,s,A,s2)
%
% Implementation of 2-mode product in 3D using loop-over-GEMMs approach
% input s not needed in our context, just for coherence with MATLAB function
% pagemtimes
  n = size(T,3);
  if(strcmp(s2,'transpose'))
    for i =1:n
      T(:,:,i) = T(:,:,i)*A.';
    end
  else
    for i =1:n
      T(:,:,i) = T(:,:,i)*A;
    end
  end
end
