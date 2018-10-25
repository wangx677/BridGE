function add_geneset(geneList,addsetName,origName,newName)

% geneList: name of genes, one list only
% addsetName: genelist name
% origName: original gene set name
% newName: assign new name to the gene set

genes = textread(geneList,'%s');
load(origName)

if size(geneset.pathwaynames,1) == 1
     pathwaynames = [geneset.pathwaynames,addsetName];
else
     pathwaynames = [geneset.pathwaynames;addsetName];
end

genenames = union(geneset.genenames,genes);

gpmatrix = zeros(length(genenames),length(pathwaynames));
for i=1:length(geneset.pathwaynames)
     tmpgenes = geneset.genenames(find(geneset.gpmatrix(:,i)==1));
     ind = find(ismember(genenames,tmpgenes));
     gpmatrix(ind,i) = 1;
end

ind = find(ismember(genenames,genes));
i = find(ismember(pathwaynames,addsetName));
gpmatrix(ind,i) = 1;

geneset.genenames = genenames';
geneset.pathwaynames = pathwaynames;
geneset.gpmatrix = gpmatrix;
geneset.entrezids = genenames';
save(newName,'geneset')

