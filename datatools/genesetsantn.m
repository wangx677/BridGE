function genesetsantn(geneListFromGeneset,geneAnnotation,outputFile)

% FUNCTION genesetsantn(geneListFromGeneset,geneAnnotation,outputFile)
%
% GENESETSANTN annotation of genes from gene sets
%
% This function is used to extract gene information from geneAnnotation,
% so that only genes from the gene sets are included.
%
% INPUTS:
%  geneListFromGeneset - a text file with each line is a gene
%  geneAnnotation - gene anotation file 
%  outputFile - output file name
%
% OUTPUTs:
%  A text file with same format as the gene anotation file
%

% read data from input files
glist = textread(geneListFromGeneset,'%s');
[chr loc1 loc2 glistall] = textread(geneAnnotation,'%s%s%s%s');

% get unique list from geneAnnotation
[C, ia, ic] = unique(glistall);
M = [chr loc1 loc2 glistall];
M = M(ia,:);

% extract genes that are included in gene sets
ind = find(ismember(M(:,4),glist)==1);
M = M(ind,:);

fwritecell(outputFile,'%s %s %s %s',M)	 
