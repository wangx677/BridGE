function msigdb2excel(symbolsFile,entrezFile)

% FUNCTION msigdb2excel(symbolsFile,entrezFile)
%
% MSIGDB2MAT convert MsigDB gene set file (.gmt) to MAT-file (MATLAB).
%
% INPUTS:
%   symbolsFile: MsigDB gene set file using gene symbols (.symbols.gmt).
%   entrezFile: MsigDB gene set file using gene entrez ids (.entrez.gmt).
%
% OUTPUTS:
%   <symbolsFile>.xlsx 

% read data files
symbolsData = readtable(symbolsFile,'FileType','text','ReadVariableNames',0,'Delimiter','\t');
entrezData = readtable(entrezFile,'FileType','text','ReadVariableNames',0,'Delimiter','\t');

for i=1:size(symbolsData,1)
     tmp = table2array(symbolsData(i,3:end));
     ind = find(cellfun(@(x)isempty(x),tmp)~=1);
     genes{i,1} = strjoin(tmp(ind),',');
     tmp = entrezData(i,ind+2);
     tmp = cellfun(@(x)num2str(x),table2cell(tmp),'uniform',0);
     entrez{i,1} = strjoin(tmp,',');
     pathwaysize(i,1) = length(ind);
end

pathway = symbolsData.Var1;
link = symbolsData.Var2;

% output file name
outputFile = strsplit(symbolsFile,'.symbols.gmt');
outputFile = sprintf('%s.xlsx',outputFile{1});

output = table(pathway,link,pathwaysize,genes,entrez);
writetable(output,outputFile);
