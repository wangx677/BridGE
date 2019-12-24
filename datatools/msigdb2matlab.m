function msigdb2matlab(symbolsFile,entrezFile)

% FUNCTION msigdb2matlab(symbolsFile,entrezFile)
% 
% MSIGDB2MAT convert MsigDB gene set file (.gmt) to MAT-file (MATLAB).
%
% INPUTS:
%   symbolsFile: MsigDB gene set file using gene symbols (.symbols.gmt).
%   entrezFile: MsigDB gene set file using gene entrez ids (.entrez.gmt).
%
% OUTPUTS:
%   <symbolsFile>.mat -This mat-file includes a strucutre array with
%      the following fields:
%      geneset.genenames - gene names
%      geneset.pathwaynames - gene set names
%      geneset.entrezids - gene entrez ids
%      geneset.gpmatrix - gene pathway binary matrix

% read data files
%symbolsData = readtable(symbolsFile,'FileType','text','ReadVariableNames',0,'Delimiter','\t');
%entrezData = readtable(entrezFile,'FileType','text','ReadVariableNames',0,'Delimiter','\t');

symbolsData = importdata(symbolsFile, ' ');
entrezData = importdata(entrezFile,' ');

% output file name
outputFile = strsplit(symbolsFile,'.symbols.gmt');
outputFile = sprintf('%s.mat',outputFile{1});

% get unique gene symbols and corresponding entrez ids
tmpsymbols = [];
tmpentrez = [];
for i=1:size(symbolsData,1)
     tmp1 = strsplit(symbolsData{i},'\t');
     tmp2 = strsplit(entrezData{i},'\t');
     pathwaynames{i} = tmp1{1};

     symbols{i} = tmp1(3:end);
     entrez{i} = tmp2(3:end);
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
