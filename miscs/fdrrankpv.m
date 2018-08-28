function fdr = fdrrankpv(pvalues,N)

% FUNCTION fdr = fdrrankpv(pvalues,N)
% FDRRANKPV compute false discovery rates (FDRs) using empirical p-values. 
% It converts pvalues to ranks and then compute FDR based on rank distributions. 
% 
% INPUTS:
% pvalues - a empirical p-value vector
% N - number of permutations
%
% OUTPUTS:
% fdr - false discovery rates
%
fdr =ones(size(pvalues));
M = unique(pvalues*N);
Mnew =  pvalues*N;
for i=1: length(M)
        ii = M(i);
        total = nnz(pvalues<=ii/N)*(ii-1)+nnz(pvalues>ii/N)*(ii);
        fdrtmp(i) =  total/N/nnz(pvalues<=ii/N);
%        fdr(Mnew==ii)=fdrtmp(i);
end

for i=1:length(M)
	ii = M(i);
	fdr(Mnew==ii)=min(fdrtmp(i:end));
end
