samplePerms = 1000;

% recessive disease model
load SNPdataAR.mat
N = length(SNPdata.pheno);

SNPdata_AR = SNPdata;
AR_SNP_protective(1,:) = hygetest(N,nnz(SNPdata_AR.pheno==0),sum(SNPdata_AR.data(find(SNPdata_AR.pheno==0),:)),sum(SNPdata_AR.data));
AR_SNP_risk(1,:) = hygetest(N,nnz(SNPdata_AR.pheno==1),sum(SNPdata_AR.data(find(SNPdata_AR.pheno==1),:)),sum(SNPdata_AR.data));

% dominant disease model
load SNPdataAD.mat
SNPdata_AD = SNPdata;
AD_SNP_protective(1,:) = hygetest(N,nnz(SNPdata_AD.pheno==0),sum(SNPdata_AD.data(find(SNPdata_AD.pheno==0),:)),sum(SNPdata_AD.data));
AD_SNP_risk(1,:) = hygetest(N,nnz(SNPdata_AD.pheno==1),sum(SNPdata_AD.data(find(SNPdata_AD.pheno==1),:)),sum(SNPdata_AD.data));

% combined disease model
combined_SNP_protective(1,:) = max(AR_SNP_protective(1,:),AD_SNP_protective(1,:));
combined_SNP_risk(1,:) = max(AR_SNP_risk(1,:),AD_SNP_risk(1,:));

% run random permutation
for R=1:samplePerms
     rng(R);
     rand_pheno = SNPdata_AR.pheno(randperm(N));
     AR_SNP_protective(R+1,:) = hygetest(N,nnz(rand_pheno==0),sum(SNPdata_AR.data(find(rand_pheno==0),:)),sum(SNPdata_AR.data));
     AR_SNP_risk(R+1,:) = hygetest(N,nnz(rand_pheno==1),sum(SNPdata_AR.data(find(rand_pheno==1),:)),sum(SNPdata_AR.data));

     AD_SNP_protective(R+1,:) = hygetest(N,nnz(rand_pheno==0),sum(SNPdata_AD.data(find(rand_pheno==0),:)),sum(SNPdata_AD.data));
     AD_SNP_risk(R+1,:) = hygetest(N,nnz(rand_pheno==1),sum(SNPdata_AD.data(find(rand_pheno==1),:)),sum(SNPdata_AD.data));

     combined_SNP_protective(R+1,:) = max(AR_SNP_protective(R+1,:),AD_SNP_protective(R+1,:));
     combined_SNP_risk(R+1,:) = max(AR_SNP_risk(R+1,:),AD_SNP_risk(R+1,:));
end

% SNP level stats
% covert to p_value
for i=1:(samplePerms+1)
     ind = setdiff(1:(samplePerms+1),i);
     AR_SNP_protective_pv(i,:) = (sum(AR_SNP_protective(ind,:)>=AR_SNP_protective(i,:))+1)/(samplePerms+1);
     AR_SNP_risk_pv(i,:) = (sum(AR_SNP_risk(ind,:)>=AR_SNP_risk(i,:))+1)/(samplePerms+1);

     AD_SNP_protective_pv(i,:) = (sum(AD_SNP_protective(ind,:)>=AD_SNP_protective(i,:))+1)/(samplePerms+1);
     AD_SNP_risk_pv(i,:) = (sum(AD_SNP_risk(ind,:)>=AD_SNP_risk(i,:))+1)/(samplePerms+1);

     combined_SNP_protective_pv(i,:) = (sum(combined_SNP_protective(ind,:)>=combined_SNP_protective(i,:))+1)/(samplePerms+1);
     combined_SNP_risk_pv(i,:) = (sum(combined_SNP_risk(ind,:)>=combined_SNP_risk(i,:))+1)/(samplePerms+1);
end

% calculate FDR
for i=1:(samplePerms+1)
     AR_SNP_protective_fdr(i,:) = mafdr(AR_SNP_protective_pv(i,:),'BH',true); 
     AR_SNP_risk_fdr(i,:) =  mafdr(AR_SNP_risk_pv(i,:),'BH',true);

     AD_SNP_protective_fdr(i,:) = mafdr(AD_SNP_protective_pv(i,:),'BH',true);
     AD_SNP_risk_fdr(i,:) = mafdr(AD_SNP_risk_pv(i,:),'BH',true);

     combined_SNP_protective_fdr(i,:) = mafdr(combined_SNP_protective_pv(i,:),'BH',true);
     combined_SNP_risk_fdr(i,:) = mafdr(combined_SNP_risk_pv(i,:),'BH',true);
end


% gene-level stats
load snpgenemapping_50kb.mat
G = length(snp2gene.genelist);
for R=0:samplePerms
     AR_gene_protective(R+1,:) = arrayfun(@(x)mean(AR_SNP_protective(R+1,find(snp2gene.sgmatrix(:,x)==1))),1:G);
     AR_gene_risk(R+1,:) = arrayfun(@(x)mean(AR_SNP_risk(R+1,find(snp2gene.sgmatrix(:,x)==1))),1:G);

     AD_gene_protective(R+1,:) = arrayfun(@(x)mean(AD_SNP_protective(R+1,find(snp2gene.sgmatrix(:,x)==1))),1:G);
     AD_gene_risk(R+1,:) = arrayfun(@(x)mean(AD_SNP_risk(R+1,find(snp2gene.sgmatrix(:,x)==1))),1:G);

     combined_gene_protective(R+1,:) = arrayfun(@(x)mean(combined_SNP_protective(R+1,find(snp2gene.sgmatrix(:,x)==1))),1:G);
     combined_gene_risk(R+1,:) = arrayfun(@(x)mean(combined_SNP_risk(R+1,find(snp2gene.sgmatrix(:,x)==1))),1:G);
end

% convert to p_value
for i=1:(samplePerms+1)
     ind = setdiff(1:(samplePerms+1),i);
     AR_gene_protective_pv(i,:) = (sum(AR_gene_protective(ind,:)>=AR_gene_protective(i,:))+1)/(samplePerms+1);
     AR_gene_risk_pv(i,:) = (sum(AR_gene_risk(ind,:)>=AR_gene_risk(i,:))+1)/(samplePerms+1);

     AD_gene_protective_pv(i,:) = (sum(AD_gene_protective(ind,:)>=AD_gene_protective(i,:))+1)/(samplePerms+1);
     AD_gene_risk_pv(i,:) = (sum(AD_gene_risk(ind,:)>=AD_gene_risk(i,:))+1)/(samplePerms+1);

     combined_gene_protective_pv(i,:) = (sum(combined_gene_protective(ind,:)>=combined_gene_protective(i,:))+1)/(samplePerms+1);
     combined_gene_risk_pv(i,:) = (sum(combined_gene_risk(ind,:)>=combined_gene_risk(i,:))+1)/(samplePerms+1);
end

% calculate FDR
for i=1:(samplePerms+1)
     AR_gene_protective_fdr(i,:) = mafdr(AR_gene_protective_pv(i,:),'BH',true);
     AR_gene_risk_fdr(i,:) = mafdr(AR_gene_risk_pv(i,:),'BH',true);

     AD_gene_protective_fdr(i,:) = mafdr(AD_gene_protective_pv(i,:),'BH',true);
     AD_gene_risk_fdr(i,:) = mafdr(AD_gene_risk_pv(i,:),'BH',true);

     combined_gene_protective_fdr(i,:) = mafdr(combined_gene_protective_pv(i,:),'BH',true);
     combined_gene_risk_fdr(i,:) = mafdr(combined_gene_risk_pv(i,:),'BH',true);
end

AR_SNP_protective = AR_SNP_protective(1,:);
AR_SNP_risk = AR_SNP_risk(1,:);
AD_SNP_protective = AD_SNP_protective(1,:);
AD_SNP_risk = AD_SNP_risk(1,:);


save('stats_single_association_and_pairwise_interaction.mat','AR*pv*','AR*fdr','AD*pv*','AD*fdr','combined*pv','combined*fdr','AR_SNP_protective','AR_SNP_risk','AD_SNP_protective','AD_SNP_risk''-v7.3')

% SNP-SNP interaction level stats
for R=0:samplePerms
     load(sprintf('ssM_mhygeSSI_alpha10.05_alpha20.05_RR_R%s.mat',num2str(R)));
     RR_ssM_protective(R+1,:) = ssM{1};
     RR_ssM_risk(R+1,:) = ssM{2};

     load(sprintf('ssM_mhygeSSI_alpha10.05_alpha20.05_DD_R%s.mat',num2str(R)));
     DD_ssM_protective(R+1,:) = ssM{1};
     DD_ssM_risk(R+1,:) = ssM{2};

     load(sprintf('ssM_mhygeSSI_alpha10.05_alpha20.05_RD_R%s.mat',num2str(R)));
     RD_ssM_protective(R+1,:) = ssM{1};
     RD_ssM_risk(R+1,:) = ssM{2};

     load(sprintf('ssM_mhygeSSI_alpha10.05_alpha20.05_combined_R%s.mat',num2str(R)));
     combined_ssM_protective(R+1,:) = ssM{1};
     combined_ssM_risk(R+1,:) = ssM{2};
end

% convert to p-value
for i=1:(samplePerms+1)
     ind = setdiff(1:(samplePerms+1),i);
     RR_ssM_protective_pv(i,:) = (sum(RR_ssM_protective(ind,:)>=RR_ssM_protective(i,:))+1)/(samplePerms+1); 
     RR_ssM_risk_pv(i,:) = (sum(RR_ssM_risk(ind,:)>=RR_ssM_risk(i,:))+1)/(samplePerms+1);

     DD_ssM_protective_pv(i,:) = (sum(DD_ssM_protective(ind,:)>=DD_ssM_protective(i,:))+1)/(samplePerms+1);
     DD_ssM_risk_pv(i,:) = (sum(DD_ssM_risk(ind,:)>=DD_ssM_risk(i,:))+1)/(samplePerms+1);

     RD_ssM_protective_pv(i,:) = (sum(RD_ssM_protective(ind,:)>=RD_ssM_protective(i,:))+1)/(samplePerms+1);
     RD_ssM_risk_pv(i,:) = (sum(RD_ssM_risk(ind,:)>=RD_ssM_risk(i,:))+1)/(samplePerms+1);

     combined_ssM_protective_pv(i,:) = (sum(combined_ssM_protective(ind,:)>=combined_ssM_protective(i,:))+1)/(samplePerms+1);
     combined_ssM_risk_pv(i,:) = (sum(combined_ssM_risk(ind,:)>=combined_ssM_risk(i,:))+1)/(samplePerms+1);
end

save('stats_single_association_and_pairwise_interaction.mat','RR*pv','RR*fdr','DD*pv','DD*fdr','RD*pv','RD*fdr','combined*pv','combined*fdr','-append','-v7.3')
