function newdata = imputesnp(data)

% FUNCTION newdata = imputesnp(data)
%
% INPUTESNP impute SNP data (0,1,2 format)
%
% This function imputes missing alleles in SNP dataset by
% replacing them with major-major.
%
% INPUTS: 
%   data - SNP data matrix (0,1,2 format, -1 is missing value) 
%     with rows are samples and columns are SNPs
%
% OUTPUTS:
%   newdata - new SNP data matrix
%

[nsample nsnp] = size(data);
for i = 1:nsnp

	snp = data(:, i);
	sum0 = sum( snp==0 );
	sum1 = sum( snp==1 );	
	sum2 = sum( snp==2 );
	[m gen] = max([sum0 sum1 sum2]);
	snp(snp==-1) = gen-1;
	data(:, i) = snp;

end
newdata = data;
