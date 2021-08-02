function computeFDR(ssmFile,samplePerms)

% Inputs
% ssmFile: SNP-SNP interaction file
% samplePerms: number of sample permutation
%
% Outputs:

for i=0:samplePerms
     load(sprintf('genstats_%s_R%s.mat',ssmFile,num2str(i)))
     % compute FDR
     % sample permutation p-value based on average interaction scores with filtering out gene-gene interactions with ranksum p>0.05
     % one tailled test
     fdr_density_protective = mafdr([bpm_pv_density{1} wpm_pv_density{1}],'BHFDR',true);
     fdr_density_risk = mafdr([bpm_pv_density{2} wpm_pv_density{2}],'BHFDR',true);     

     % two tailed test
     fdr_tmp = mafdr([bpm_pv_density{1} wpm_pv_density{1} bpm_pv_density{2} wpm_pv_density{2}],'BHFDR',true);
     fdr_density_protective_both = fdr_tmp(1:length(fdr_tmp)/2);
     fdr_density_risk_both = fdr_tmp(length(fdr_tmp)/2+1:end);

     % compute FDR
     % sample permutation p-value based on average interaction scores and ranksum score
     % one tailled test
     fdr_density_and_ranksum_protective = mafdr([bpm_pv_density_and_ranksum{1} wpm_pv_density_and_ranksum{1}],'BHFDR',true);
     fdr_density_and_ranksum_risk = mafdr([bpm_pv_density_and_ranksum{2} wpm_pv_density_and_ranksum{2}],'BHFDR',true);

     % two tailed test
     fdr_tmp = mafdr([bpm_pv_density_and_ranksum{1} wpm_pv_density_and_ranksum{1} bpm_pv_density_and_ranksum{2} wpm_pv_density_and_ranksum{2}],'BHFDR',true);
     fdr_density_and_ranksum_protective_both = fdr_tmp(1:length(fdr_tmp)/2);
     fdr_density_and_ranksum_risk_both = fdr_tmp(length(fdr_tmp)/2+1:end);

     save(sprintf('genstats_%s_R%s.mat',ssmFile,num2str(i)),'fdr*','-append')     
end
