function [phenonew,famidtmp,indidtmp]  = withinclassrand(datafile,plinkCluster2,randseed)

% FUNCTION [phenonew,famidtmp,indidtmp]  = withinclassrand(datafile,plinkCluster2,randseed)
% 
% WITHINCLASSRAND generates random pheno lables based on matched cases and controls.
%
% INPUTS:
%    datafile - MATLAB data file, can be SNPdataAR.mat, SNPdataAD.mat, or gwas_data_final.mat
%    plinkCluster2 - output from plink size-2 clustering after removing individuals that
%         are not paired with others.
%    randseed - random seed
% OUTPUTS:
%    phenonew - randomized phenotype vector
%

load(datafile)
pheno_orig = SNPdata.pheno;
[famidtmp indidtmp clustertmp] = textread(plinkCluster2,'%s%s%s');
ind2keep = find(ismember(famidtmp,SNPdata.fid)==1);

famidtmp = famidtmp(ind2keep);
indidtmp = indidtmp(ind2keep);
clustertmp = clustertmp(ind2keep);

clustertmpuniq = unique(clustertmp); 

for i=1:length(clustertmpuniq)
	ind = find(ismember(clustertmp,clustertmpuniq(i))==1);
	if length(ind)==2
		clusterpair(i,:) = [famidtmp(ind(1)) indidtmp(ind(1)) clustertmp(ind(1)) famidtmp(ind(2)) indidtmp(ind(2)) clustertmp(ind(2))];
	else
		sprintf('cluster pair error')
		return
	end
end

rand('seed',randseed);

randidx = randperm(size(clusterpair,1));

n = floor(size(clusterpair,1)/2);

f1 = clusterpair(:,1);
p1 = clusterpair(:,2);
f2 = clusterpair(:,4);
p2 = clusterpair(:,5);

phenonew = SNPdata.pheno;

ind = find(ismember(SNPdata.fid,f1(randidx(1:n)))==1 & ismember(SNPdata.pid,p1(randidx(1:n)))==1);
phenonew(ind)  = 1 - SNPdata.pheno(ind);

ind = find(ismember(SNPdata.fid,f2(randidx(1:n)))==1 & ismember(SNPdata.pid,p2(randidx(1:n)))==1);
phenonew(ind)  = 1 - SNPdata.pheno(ind); 
