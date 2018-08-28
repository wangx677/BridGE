function gmt2matlab(symbolsFile,entrezFile)

% FUNCTION gmt2mat(symbolsFile,entrezFile)
% 
% MSIGDB2MAT convert MsigDB gene set file (.gmt) to MAT-file (MATLAB).
%
% INPUTS:
%   symbolsFile: gmt file using gene symbols
%
% OUTPUTS:
%   <symbolsFile>.mat -This mat-file includes a strucutre array with
%      the following fields:
%      geneset.genenames - gene names
%      geneset.pathwaynames - gene set names
%      geneset.gpmatrix - gene pathway binary matrix

% rename files to csv
symbolsFile = strsplit(symbolsFile,'.gmt');
symbolsFile = symbolsFile{1};

% output file name
outputFile = strsplit(symbolsFile,'.symbols');
outputFile = sprintf('%s.mat',outputFile{1});

% read data files
symbolsData = importdata(sprintf('%s.gmt',symbolsFile));


% get unique gene symbols and corresponding entrez ids
tmpsymbols = [];
for i=1:length(symbolsData)
     tmp = strsplit(symbolsData{i},'\t');
     pathwaynames{i} = tmp{2};
     symbols{i} = tmp(3:end);
     tmpsymbols = [tmpsymbols symbols{i}];
end

[genenames ia ic] = unique(tmpsymbols);

% binary matrix for gene and pathway association
gpmatrix = zeros(length(genenames),length(pathwaynames));

for i=1:length(pathwaynames)
     ind = find(ismember(genenames,symbols{i}));
     gpmatrix(ind,i) = 1;
end

% save data as structure array
geneset.genenames = genenames;
geneset.pathwaynames = pathwaynames;
geneset.gpmatrix = gpmatrix;
save(sprintf('%s',outputFile),'geneset','-v7.3');
