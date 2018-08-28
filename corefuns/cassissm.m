function cassissm(cassiFile)

% FUNCTION cassissm(cassiFile)
%
% CASSISSM converts SNP-SNP interactions that are computed using CASSI
%
% INPUTS:
%   cassiFile - outputfile from CASSI
% 
% OUTPUTS:
%   A mat-file include a matrix vector (size of 2) ssM. ssM{1} is for protective interactions and
%   ssM{2} is for risk interactions.
%

load SNPdataAD.mat
rsid = SNPdata.rsid;

C = strsplit(cassiFile,'.');
p = length(SNPdata.rsid);
q = p;

if ismember(C(end),'lr')==1
	% logistic regression (cassi)
     system(sprintf('mv %s %s.csv',cassiFile,cassiFile));
     data = readtable(sprintf('%s.csv',cassiFile));
     x = data.SNP1;
     y = data.SNP2;
     w = data.LR_LOG_OR;
     z = data.LR_P;
     snp1 = data.ID1;
     snp2 = data.ID2;
	clear data
end

system(sprintf('mv %s.csv %s',cassiFile,cassiFile));

n = [x;y];
snps = [snp1;snp2]; 
clear snp1 snp2
[un in]=unique(n);
clear n

z = -log10(z);
cassi = sparse(x,y,z,p,q);
if exist('log_or')~=1
        log_or = sparse(x,y,w,p,q);
end

clear x y w p q

ssM{1} = cassi.*(log_or<0);
ssM{2} = cassi.*(log_or>0);

clear cassi log_or

for tt=1:2
	ssM{tt} = triu(ssM{tt},1)+triu(ssM{tt},1)';
end


if isequal(snps(in),rsid(un))~=1
	[tmp1 rs tmp3 tmp4 tmp5 tmp6]=textread(sprintf('%s/gwas_data_final.bim',inputdir),'%s%s%s%s%s%s');
	if isequal(snps(in),rs(un))==1;
		[tmp ind]=ismember(rs,rsid);
		for tt=1:2
			ssM{tt} = ssM{tt}(ind,ind);
		end
	else
		for tt=1:2
			ssM{tt} = [];
		end
	end
end
save(sprintf('ssM_%s_%s.%s.mat',C{3},C{1},C{2}),'ssM','-v7.3')
