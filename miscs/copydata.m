function copydata(old,new)

% copy data from old BridGE to new BridGE
% old = 'Diabetes1_EGA08';
% new = 'project_Diabetes1_EGA08';

% gwas_data_final.*

system(sprintf('cp /project/chadm/wwang/BPM2015/%s/gwas_data_final.* /project/csbio/wwang/BridGE/%s/',old,new));

cd(sprintf('/project/csbio/wwang/BridGE/%s/',new))

load(sprintf('/project/chadm/wwang/BPM2015/%s/gwas_data_final.mat',old))

SNPdata.data = trans_data;
SNPdata.rsid = snpRSNums;
SNPdata.chr = chromoID;
SNPdata.loc = chromoLoc'
SNPdata.pheno = pheno;
SNPdata.pid = patient_ids;
SNPdata.fid = family_ids;
SNPdata.gender = gender;

save('gwas_data_final.mat','SNPdata')

bindataar('gwas_data_final.mat');
bindataad('gwas_data_final.mat');

clearvars -except old new
% cluster2 file
system(sprintf('cp /project/chadm/wwang/BPM2015/%s/PlinkFile.cluster2 /project/csbio/wwang/BridGE/%s/PlinkFile.cluster2',old,new));

% snp_pathway_min10_max300.mat
load(sprintf('/project/chadm/wwang/BPM2015/%s/snp_pathway_min10_max300.mat',old))
snpset.spmatrix = snp_pathway;
snpset.pathwaynames = pathwayNames;
snpset.snplist = snpList;
snpset.pathsize = sum(snp_pathway)';
snpset.genesets = 'c2.all.v3.0.symbols_CP.mat';
save('snp_pathway_min10_max300.mat','snpset');
clearvars -except old new

% snpgenemapping_50kb.mat
load(sprintf('/project/chadm/wwang/BPM2015/%s/snpgenemapping_50kb.mat',old))
snp2gene.snplist = snpList;
snp2gene.genelist = geneList;
snp2gene.sgmatrix = snp_gene_matrix;
save('snpgenemapping_50kb.mat','snp2gene');
clearvars -except old new

% BPMind.mat
bpmind('snp_pathway_min10_max300.mat')


