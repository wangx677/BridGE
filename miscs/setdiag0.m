function M = setdiag0(M);

%SETDIAG0 set the diagonal of a symmetric matrix to be all 0s
%
% INPUT:
%   M: a symmetric matrix

if (isequal(triu(M,1),tril(M,-1)')==1)
     M = triu(M,1)+tril(M,-1);
else
     error('Input matrix is not symmetric!')
end
