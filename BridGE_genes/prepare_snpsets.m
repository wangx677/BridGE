function prepare_snpsets(gene_annotation,mappingDistance)

% get snp to gene mapping
if mappingDistance<1000
     snpGeneMappingFile = sprintf('snpgenemapping_%sbp.mat',num2str(mappingDistance));
else
     snpGeneMappingFile = sprintf('snpgenemapping_%skb.mat',num2str(mappingDistance/1000));
end
bimfile='gwas_data_final.bim';
mapsnp2gene(bimfile,gene_annotation,mappingDistance,'matrix',snpGeneMappingFile)

% create gene set
tmp = readtable('gene_annotation','filetype','text');
geneset.genenames = tmp.Var4;
geneset.pathwaynames = tmp.Var4;
geneset.gpmatrix = eye(length(tmp.Var4),length(tmp.Var4));
save('GOI_geneset','geneset')

% get snp to pbody gene mapping
snppathway('SNPdataAR.mat', snpGeneMappingFile,'GOI_geneset.mat',5,500)

% get BPMind file (snp set for each gene)
bpmind('snp_pathway_min5_max500.mat')
