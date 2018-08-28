function genesetglist(genesets)

% FUNCTION genesetglist(genesets)
%
% GENESETGLIST print the list of genes from gene set
%
% INPUTS:
%   genesets - gene set file name. It is a mat-file includes a strucutre array with
%      the following fields:
%      geneset.genenames - gene names
%      geneset.pathwaynames - gene set names
%      geneset.entrezids - gene entrez ids
%      geneset.gpmatrix - gene pathway binary matrix
% OUTPUTS:
%   'geneList_from_genesets' - a text file includes the list of genes from gene sets  

load(genesets)
fwritecell('geneList_from_genesets','%s',geneset.genenames)
