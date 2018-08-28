function msigdb2matlab(symbolsFile,entrezFile)

% FUNCTION msigdb2mat(symbolsFile,entrezFile)
% 
% MSIGDB2MAT convert MsigDB gene set file (.gmt) to MAT-file (MATLAB).
%
% INPUTS:
%   symbolsFile: MsigDB gene set file using gene symbols (.symbols.gmt).
%   entrezFile: MsigDB gene set file using gene entrez ids (.symbols.gmt).
%
% OUTPUTS:
%   <symbolsFile>.mat -This mat-file includes a strucutre array with
%      the following fields:
%      geneset.genenames - gene names
%      geneset.pathwaynames - gene set names
%      geneset.entrezids - gene entrez ids
%      geneset.gpmatrix - gene pathway binary matrix

% rename files to csv
symbolsFile = strsplit(symbolsFile,'.gmt');
symbolsFile = symbolsFile{1};
system(sprintf('cp %s.gmt %s.csv',symbolsFile,symbolsFile));

entrezFile = strsplit(entrezFile,'.gmt');
entrezFile = entrezFile{1};
system(sprintf('cp %s.gmt %s.csv',entrezFile,entrezFile));

% output file name
outputFile = strsplit(symbolsFile,'.symbols');
outputFile = sprintf('%s.mat',outputFile{1});

% read data files
symbolsData = readtable(sprintf('%s.csv',symbolsFile),'headerlines',0);
delete(sprintf('%s.csv',symbolsFile));

entrezData = readtable(sprintf('%s.csv',entrezFile),'headerlines',0);
delete(sprintf('%s.csv',entrezFile));

% only need Var7 from the data
symbolsData = symbolsData.Var7;
entrezData = entrezData.Var7;

% get unique gene symbols and corresponding entrez ids
tmpsymbols = [];
tmpentrez = [];
for i=1:length(symbolsData)
     tmp1 = strsplit(symbolsData{i});
     tmp2 = strsplit(entrezData{i});
     pathwaynames{i} = tmp1{1};

     symbols{i} = tmp1(2:end);
     entrez{i} = tmp2(2:end);
     tmpsymbols = [tmpsymbols symbols{i}];
     tmpentrez = [tmpentrez entrez{i}];
end

[genenames ia ic] = unique(tmpsymbols);
entrezids = tmpentrez(ia);

% binary matrix for gene and pathway association
gpmatrix = zeros(length(genenames),length(pathwaynames));

for i=1:length(pathwaynames)
     ind = find(ismember(genenames,symbols{i}));
     gpmatrix(ind,i) = 1;
end

% save data as structure array
geneset.genenames = genenames;
geneset.pathwaynames = pathwaynames;
geneset.entrezids = entrezids;
geneset.gpmatrix = gpmatrix;
save(sprintf('%s',outputFile),'geneset','-v7.3');
