cd /project/csbio/wwang/BridGE_nobin/project_PD_Simon
load genstats_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat
pv_perm = [bpm_local_pv{1} bpm_local_pv{2}];
load genstats_zscore_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0_snpPerm1000.mat
pv_zscore = [bpm_pv_snp{1} bpm_pv_snp{2}];
ind = find(pv_zscore~=1 & pv_perm~=1);

subplot(2,2,1)
scatter(-log10(pv_perm(ind)),-log10(pv_zscore(ind)),'.')
xlabel('permutation p-value (-log10)')
ylabel('zscore p-value (-log10)')
title(sprintf('PD Simon1 corr=%.2f',corr(-log10(pv_perm(ind))',-log10(pv_zscore(ind))')))
setfig

cd /project/csbio/wwang/BridGE_nobin/project_PD_Simon_test_updated_ssM/
load genstats_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat
pv_perm = [bpm_local_pv{1} bpm_local_pv{2}];
load genstats_zscore_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0_snpPerm1000.mat
pv_zscore = [bpm_pv_snp{1} bpm_pv_snp{2}];
ind = find(pv_zscore~=1 & pv_perm~=1);
 
subplot(2,2,2)
scatter(-log10(pv_perm(ind)),-log10(pv_zscore(ind)),'.')
xlabel('permutation p-value (-log10)')
ylabel('zscore p-value (-log10)')
title(sprintf('PD Simon2 corr=%.2f',corr(-log10(pv_perm(ind))',-log10(pv_zscore(ind))')))

cd /project/csbio/wwang/BridGE_nobin/project_PD_NGRC
load genstats_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat
pv_perm = [bpm_local_pv{1} bpm_local_pv{2}];
load genstats_zscore_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0_snpPerm1000.mat
pv_zscore = [bpm_pv_snp{1} bpm_pv_snp{2}];
ind = find(pv_zscore~=1 & pv_perm~=1);
 
setfig
subplot(2,2,3)
scatter(-log10(pv_perm(ind)),-log10(pv_zscore(ind)),'.')
xlabel('permutation p-value (-log10)')
ylabel('zscore p-value (-log10)')
title(sprintf('PD NGRC1 corr=%.2f',corr(-log10(pv_perm(ind))',-log10(pv_zscore(ind))')))
setfig

cd /project/csbio/wwang/BridGE_nobin/project_PD_NGRC_test_updated_ssM/
load genstats_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat
pv_perm = [bpm_local_pv{1} bpm_local_pv{2}];
load genstats_zscore_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0_snpPerm1000.mat
pv_zscore = [bpm_pv_snp{1} bpm_pv_snp{2}];
subplot(2,2,4)
ind = find(pv_zscore~=1 & pv_perm~=1);
 
scatter(-log10(pv_perm(ind)),-log10(pv_zscore(ind)),'.')
xlabel('permutation p-value (-log10)')
ylabel('zscore p-value (-log10)')
title(sprintf('PD NGRC2 corr=%.2f',corr(-log10(pv_perm(ind))',-log10(pv_zscore(ind))')))
setfig

cd /project/csbio/wwang/BridGE_nobin
saveas(gcf,'scatter_comparison_pv_perm_vs_zscore.jpg')
saveas(gcf,'scatter_comparison_pv_perm_vs_zscore.pdf')

