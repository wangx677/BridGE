function update_pathway(genes,snpgenemapping,snppathwayfile,genesetname)

% this script is used to update the snppathwayfile when testing for interaction
% between a new geneset and the existing pathways

load(genes)
load(snpgenemapping)

% get snp index that can be mapped to genes
gene_ind = find(ismember(snp2gene.genelist,genes));
gene_ind = find(sum(snp2gene.sgmatrix(:,gene_ind),2)~=0);

load(snppathwayfile)
load(sprintf('%s/refdata/%s',getenv('BRIDGEPATH'),snpset.genesets))

snpset.spmatrix(:,size(snpset.spmatrix,2)+1) = zeros(size(snpset.spmatrix,1),1);
snpset.spmatrix(gene_ind,end) = 1;
snpset.pathsize(end+1) = length(gene_ind);
snpset.pathwaynames{end+1} = genesetname;


gene_ind = find(ismember(geneset.genenames,genes));
geneset.gpmatrix(:,size(geneset.gpmatrix,2)+1) = zeros(size(geneset.gpmatrix,1),1);
geneset.gpmatrix(gene_ind,end) = 1;
geneset.pathwaynames{end+1} = genesetname;
 
newgeneset = strsplit(snpset.genesets,'.mat');
newgeneset = sprintf('%s.%s.mat',newgeneset{1},genesetname);
save(sprintf('%s/refdata/%s',getenv('BRIDGEPATH'),newgeneset),'geneset')

snpset.genesets = newgeneset;

newsnpset = strsplit(snppathwayfile,'.mat');
newsnpset = sprintf('%s_%s.mat',newsnpset{1},genesetname);
save(newsnpset,'snpset');
