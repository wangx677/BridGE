function outputFile = snppathway(dataFile,snpGeneMappingFile,genesets,minPath,maxPath)

% FUNCTION outputFile = snppathway(dataFile,snpGeneMappingFile,genesets,minPath,maxPath)
%
% SNPPATHWAY construct snp pathway association based on given snps, snp-gene mapping, 
%  and gene-pathway mapping.
%
% INPUTS:
%   dataFile - a mat-file includes a structure array SNPdata with fields of data, rsid, chr, 
%     loc, pheno, fid, pid, and gender;
%   snpGeneMappingFile - a mat-file includes a structure array snp2gene with fields of sgmatrix, 
%     snplist, and genelist
%   genesets - a mat-file includes a structure array geneset with fields of genenames, 
%     pathwaynames, entrezids (optional), and GenePathwayMatrix
%   minPath - minimum pathway size
%   maxPath - maximum pathway size
%
% OUTPUTS:
%   snp_pathway_min<minPath>_max<maxPath>.mat, a mat-file includes a structure array snpset with the 
%     following fields:
%     spmatrix - a binary matrix showing SNP-pathway relationships. Rows are SNPs and columns are pathways
%     pathwaynames - pathway names
%     snplist - SNP ids
%     pathsize - pathway sizes
%

% load file
% snpGeneMappingFile includes snp2gene.snplist, snp2gene.genelist, snp2gene.sgmatrix
load(snpGeneMappingFile); 

% dataFile includes SNPdata.rsid, SNPdata.data, SNPdata.chr, SNPdata.loc, SNPdata.pheno, SNPdata.fid, SNPdata.pid, SNPdata.gender
load(dataFile);

% genesets includes geneset.genenames, geneset.pathwaynames, geneset.entrezids(optional), geneset.gpmatrix
load(genesets); 

% only keep genes we have in the data
ind = ismember(geneset.genenames,snp2gene.genelist);
geneset.gpmatrix  =  geneset.gpmatrix(ind,:);
geneset.genenames = geneset.genenames(ind);
if isfield(geneset,'entrezids')
     geneset.entrezids = geneset.entrezids(ind);
end

% remove pathways that with size less than minPath or greater than maxPath
% ind = find(sum(geneset.gpmatrix)>=minPath & sum(geneset.gpmatrix)<=maxPath);
% geneset.gpmatrix = geneset.gpmatrix(:,ind);
% geneset.pathwaynames = geneset.pathwaynames(ind);

% remove genes that can't be mapped to any pathways
% ind = find(sum(geneset.gpmatrix,2)==0);
% geneset.genenames(ind)=[];
% geneset.gpmatrix(ind,:)=[];
% if isfield(geneset,'entrezids')
%      geneset.entrezids(ind) = [];
% end

% reorder geneset.gpmatrix and snp2gene.sgmatrix to make genes in two matrices have same order
[t inda indb] = intersect(geneset.genenames,snp2gene.genelist);

geneset.genenames = geneset.genenames(inda);
if isfield(geneset,'entrezids')
     geneset.entrezids = geneset.entrezids(inda);
end
geneset.gpmatrix = geneset.gpmatrix(inda,:);

snp2gene.genelist = snp2gene.genelist(indb);
snp2gene.sgmatrix = snp2gene.sgmatrix(:,indb);

% reorder snp2gene.sgmatrix to make snps in the matrix have same order as SNPdata.rsid in the data

[t ind] = ismember(SNPdata.rsid,snp2gene.snplist);
snp2gene.snplist(ind==0)=[];
snp2gene.sgmatrix(ind==0,:)=[];
[t ind] = ismember(SNPdata.rsid,snp2gene.snplist);
snp2gene.snplist = snp2gene.snplist(ind);
snp2gene.sgmatrix = snp2gene.sgmatrix(ind,:);


% get snp-pathway matrix

spmatrix = (snp2gene.sgmatrix*1)*geneset.gpmatrix;

spmatrix = spmatrix>0;

ind = find(sum(spmatrix)>=minPath & sum(spmatrix)<=maxPath);
spmatrix = spmatrix(:,ind);
geneset.pathwaynames = geneset.pathwaynames(ind);

pathsize = sum(spmatrix);

snpset.spmatrix = spmatrix;
snpset.pathwaynames = geneset.pathwaynames;
snpset.snplist = snp2gene.snplist;
snpset.pathsize = pathsize;
genesets = strsplit(genesets,'/');
genesets = genesets{end};
snpset.genesets = genesets;

outputFile=sprintf('snp_pathway_min%d_max%d.mat',minPath,maxPath);
save(outputFile, 'snpset','-v7.3')
end
