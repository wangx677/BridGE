outputFile = 'c2.all.v6.0.symbols_hap';
data = readtable('Gene_level_statistics_only_library_genes.xlsx');
pathwaynames = unique(data.Complex_PW);

% get pathway genes
tmpsymbols = [];
for i=1:length(pathwaynames)
     symbols{i} = data.Genes_in_Complex_PW(find(ismember(data.Complex_PW,pathwaynames{i})));
     tmpsymbols = [tmpsymbols symbols{i}'];
end

% get unique gene names
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
save(sprintf('%s.mat',outputFile),'geneset','-v7.3');

