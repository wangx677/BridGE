function summarize_individual_result_add_ld(GI,R)
% get gene-gene interactions that have FDR<=0.25 or pv<=0.05

% GI: Genetic interaction type: 'hygeSSI' or 'mhygeSSI';
% R: 0 is for real label, 1, 2, 3, ... for random labels


model={'RR','DD','RD','combined'};

load BPMind.mat
BPMsize = reshape(BPM.size,length(BPM.size),1);
WPMsize = reshape(WPM.size,length(WPM.size),1);

% get gene names for between-gene interaction and within-gene interaction
BPM_gene1 = WPM.pathway(BPM.path1idx);
BPM_gene2 = WPM.pathway(BPM.path2idx);
WPM_gene1 = WPM.pathway;
WPM_gene2 = WPM.pathway;

k=1;
% collect results from 4 different disease models: RR, DD, RD, combined
for i=1:length(model)
     load(sprintf('genstats_ssM_%s_alpha10.05_alpha20.05_%s_R%s.mat',GI,model{i},num2str(R)))
     pv_protective(:,k) = [bpm_pv_density_and_ranksum{1}';wpm_pv_density_and_ranksum{1}'];
     pv_risk(:,k) = [bpm_pv_density_and_ranksum{2}';wpm_pv_density_and_ranksum{2}'];

     ranksum_local_protective(:,k) = [bpm_local{1}';wpm_local{1}'];
     ranksum_local_risk(:,k) = [bpm_local{2}';wpm_local{2}'];
     ranksum_local_protective_ld(:,k) = [bpm_local_ld{1}';wpm_local_ld{1}'];
     ranksum_local_risk_ld(:,k) = [bpm_local_ld{2}';wpm_local_ld{2}'];

     
     density_protective(:,k) = [BPM_density{1}';WPM_density{1}'];
     density_risk(:,k) = [BPM_density{2}';WPM_density{2}'];
     density_protective_ld(:,k) = [BPM_density_ld{1}';WPM_density_ld{1}'];
     density_risk_ld(:,k) = [BPM_density_ld{2}';WPM_density_ld{2}'];

     % FDR based on one-tailed test
     fdr_protective(:,k) = fdr_density_and_ranksum_protective';
     fdr_risk(:,k) = fdr_density_and_ranksum_risk';

     % FDR based on two-tailed test
     fdr_protective_both(:,k) = fdr_density_and_ranksum_protective_both';
     fdr_risk_both(:,k) = fdr_density_and_ranksum_risk_both';

     k = k+1;
end

% for each interaction, find the disease model that has best score
% protective interactions (one-tailed test)
% ind1 = find(min(BPM_fdr_protective,[],2)<=0.25 | min(BPM_pv_protective,[],2)<=0.05);
% ind2 = find(min(WPM_fdr_protective,[],2)<=0.25 | min(WPM_pv_protective,[],2)<=0.05)
gene1 = [BPM_gene1;WPM_gene1];
gene2 = [BPM_gene2;WPM_gene2];
[pv idx] = min(pv_protective,[],2);

for i=1:length(idx)
     fdr(i,1) = fdr_protective(i,idx(i));
     fdr_both(i,1) = fdr_protective_both(i,idx(i));
     density(i,1) = density_protective(i,idx(i));
     density_ld(i,1) = density_protective_ld(i,idx(i));
     ranksum_local(i,1) = ranksum_local_protective(i,idx(i));
     ranksum_local_ld(i,1) = ranksum_local_protective_ld(i,idx(i));
     dmodel(i,1) = model(idx(i));
end

msize = [BPMsize;WPMsize];

result_protective_all = table(gene1,gene2,dmodel,msize,density,density_ld,ranksum_local,ranksum_local_ld,pv,fdr,fdr_both);
[tmp idx] = sort(result_protective_all.pv,'ascend');
result_protective_all = result_protective_all(idx,:);

ind = find(result_protective_all.pv<=0.05);
result_protective_marginal_sig = result_protective_all(ind,:);

clear fdr density dmodel fdr_both density_ld ranksum_local ranksum_local_ld

% risk interactions (one-tailed test)
% ind1 = find(min(BPM_fdr_risk,[],2)<=0.25 | min(BPM_pv_risk,[],2)<=0.05);
% ind2 = find(min(WPM_fdr_risk,[],2)<=0.25 | min(WPM_pv_risk,[],2)<=0.05)
gene1 = [BPM_gene1;WPM_gene1];
gene2 = [BPM_gene2;WPM_gene2];
[pv idx] = min(pv_risk,[],2);

for i=1:length(idx)
     fdr(i,1) = fdr_risk(i,idx(i));
     fdr_both(i,1) = fdr_risk_both(i,idx(i));
     density(i,1) = density_risk(i,idx(i));
     density_ld(i,1) = density_risk_ld(i,idx(i));
     ranksum_local(i,1) = ranksum_local_risk(i,idx(i));
     ranksum_local_ld(i,1) = ranksum_local_risk_ld(i,idx(i));
     dmodel(i,1) = model(idx(i));
end

result_risk_all = table(gene1,gene2,dmodel,msize,density,density_ld,ranksum_local,ranksum_local_ld,pv,fdr,fdr_both);
[tmp idx] = sort(result_risk_all.pv,'ascend');
result_risk_all = result_risk_all(idx,:);

ind = find(result_risk_all.pv<=0.05);
result_risk_marginal_sig = result_risk_all(ind,:);
clear fdr density dmodel fdr density dmodel

save(sprintf('summary_results_%s_R%s.mat',GI,num2str(R)),'result_risk*','result_protective*')
