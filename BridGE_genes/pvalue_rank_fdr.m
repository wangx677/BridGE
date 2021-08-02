function [fdr] = pvalue_rank_fdr(pvalue,nn)

fdr =ones(size(pvalue));
mm = unique(pvalue*nn);
mm_new =  pvalue*nn;
for i=1: length(mm)
        ii = mm(i);
        total = nnz(pvalue<=ii/nn)*(ii-1)+nnz(pvalue>ii/nn)*(ii);
        fdrtmp(i) =  total/nn/nnz(pvalue<=ii/nn);
%        fdr(mm_new==ii)=fdrtmp(i);
end

for i=1:length(mm)
	ii = mm(i);
	fdr(mm_new==ii)=min(fdrtmp(i:end));
end
