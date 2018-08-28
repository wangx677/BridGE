function fdr = rankfdr(pvalue,nperms)

% FUNCTION fdr = rankfdr(pvalue,nperms)
% 
% RANKFDR estimates false discovery rate based on empirical p-values. 
% For each unique p from pvalue, we can calculate the number of pvalue<=p (n1).
% Since pvalue is derived from permutations, we will be able to know  the 
% number of PVALUE<=p (n2) in the random runs by the given pvalue.
% 
% For example, if the size of pvalue is N, and there are n1 p<=1/nperms in pvalue.
% then in total there will be (N-n1) p'<=1/nperms from permutations, so n2 = N-n1.
% the fdr of p=1/nperms can be calculated as (n2/nperms)/n1.
%
% This idea can be generalized to any p=j/nperms. 
% n1 = nnz(pvalue<=j/nperms);
% n2 = nnz(pvalue<=j/nperms)*(j-1)+nnz(pvalue>j/nperms)*(j);
% fdr = (n2/nperms)/n1;
%
% INPUTS:
%   pvalue - a vector of p-values
%   nperms - number of permutations used to derive the p-values
%
% OUTPUTS:
%  fdr - false discovery rate

fdr =ones(size(pvalue));
mm = unique(pvalue*nperms);
mm_new =  pvalue*nperms;
for i=1: length(mm)
        j = mm(i);
        total = nnz(pvalue<=j/nperms)*(j-1)+nnz(pvalue>j/nperms)*(j);
        fdrtmp(i) =  total/nperms/nnz(pvalue<=j/nperms);
%        fdr(mm_new==j)=fdrtmp(i);
end

% if the fdr of a p-value is greater than a larger p-value
% update this fdr to be same as the larger p-value's fdr
for i=1:length(mm)
	j = mm(i);
	fdr(mm_new==j)=min(fdrtmp(i:end));
end
