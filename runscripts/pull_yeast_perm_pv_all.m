phenos = {'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};
models = {'DD','RD','RR','combined'};

% same pathway standard was used across data

for i=1:length(phenos)
     projectdir=sprintf('%s/project_yeast_%s_complex_t25_b50_mhygeSSI',getenv('BRIDGEPATH'),phenos{i});
     load(sprintf('%s/BPMind.mat',projectdir));
     load(sprintf('%s/snp_pathway_min5_max300.mat',projectdir));
     path1 = snpset.pathwaynames(BPM.path1idx);
     path2 = snpset.pathwaynames(BPM.path2idx);
     path1 = [reshape(path1,length(path1),1);reshape(path1,length(path1),1)];
     path2 = [reshape(path2,length(path2),1);reshape(path2,length(path2),1)];

     path = snpset.pathwaynames;
     path = [reshape(path,length(path),1);reshape(path,length(path),1)];

    for j=1:length(models)
          load(sprintf('%s/results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',projectdir,models{j}));
          BPM_tmp(:,j) = bpm_pv';
          WPM_tmp(:,j) = wpm_pv';
          PATH_tmp(:,j) = path_pv';
     end
     BPM_pv(:,i) = min(BPM_tmp,[],2);
     WPM_pv(:,i) = min(WPM_tmp,[],2);
     PATH_pv(:,i) = min(PATH_tmp,[],2);
end
bpm_eff = [repmat({'protective'},size(BPM_pv,1)/2,1); repmat({'risk'},size(BPM_pv,1)/2,1)];
wpm_eff = [repmat({'protective'},size(WPM_pv,1)/2,1); repmat({'risk'},size(WPM_pv,1)/2,1)];

output_bpm_t25 = array2table(BPM_pv);
output_bpm_t25.Properties.VariableNames = phenos;
output_bpm_t25 = [table(path1,path2,bpm_eff) output_bpm_t25];
[tmp ind] = sort(sum(BPM_pv<=0.05,2),'descend');
output_bpm_t25 = output_bpm_t25(ind,:);
output_bpm_t25 = output_bpm_t25(1:nnz(tmp>=2),:);

output_wpm_t25 = array2table(WPM_pv);
output_wpm_t25.Properties.VariableNames = phenos;
output_wpm_t25 = [table(path,wpm_eff) output_wpm_t25];
[tmp ind] = sort(sum(WPM_pv<=0.05,2),'descend');
output_wpm_t25 = output_wpm_t25(ind,:);
output_wpm_t25 = output_wpm_t25(1:nnz(tmp>=2),:);

writetable(output_bpm_t25,sprintf('%s/results_collection/yeast/yeast_across_pheno_pvalues_t25.xls',getenv('BRIDGEPATH')),'Sheet',1)
writetable(output_wpm_t25,sprintf('%s/results_collection/yeast/yeast_across_pheno_pvalues_t25.xls',getenv('BRIDGEPATH')),'Sheet',2)
clear BPM_tmp WPM_tmp PATH_tmp BPM_pv WPM_pv PATH_pv
                
for i=1:length(phenos)
     projectdir=sprintf('%s/project_yeast_%s_complex_t50_b25_mhygeSSI',getenv('BRIDGEPATH'),phenos{i});
     load(sprintf('%s/BPMind.mat',projectdir));
     load(sprintf('%s/snp_pathway_min5_max300.mat',projectdir));
     path1 = snpset.pathwaynames(BPM.path1idx);
     path2 = snpset.pathwaynames(BPM.path2idx);
     path1 = [reshape(path1,length(path1),1);reshape(path1,length(path1),1)];
     path2 = [reshape(path2,length(path2),1);reshape(path2,length(path2),1)];

     path = snpset.pathwaynames;
     path = [reshape(path,length(path),1);reshape(path,length(path),1)];

    for j=1:length(models)
          load(sprintf('%s/results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',projectdir,models{j}));
          BPM_tmp(:,j) = bpm_pv';
          WPM_tmp(:,j) = wpm_pv';
          PATH_tmp(:,j) = path_pv';
     end
     BPM_pv(:,i) = min(BPM_tmp,[],2);
     WPM_pv(:,i) = min(WPM_tmp,[],2);
     PATH_pv(:,i) = min(PATH_tmp,[],2);
end
bpm_eff = [repmat({'protective'},size(BPM_pv,1)/2,1); repmat({'risk'},size(BPM_pv,1)/2,1)];
wpm_eff = [repmat({'protective'},size(WPM_pv,1)/2,1); repmat({'risk'},size(WPM_pv,1)/2,1)];

output_bpm_t50 = array2table(BPM_pv);
output_bpm_t50.Properties.VariableNames = phenos;
output_bpm_t50 = [table(path1,path2,bpm_eff) output_bpm_t50];
[tmp ind] = sort(sum(BPM_pv<=0.05,2),'descend');
output_bpm_t50 = output_bpm_t50(ind,:);
output_bpm_t50 = output_bpm_t50(1:nnz(tmp>=2),:);

output_wpm_t50 = array2table(WPM_pv);
output_wpm_t50.Properties.VariableNames = phenos;
output_wpm_t50 = [table(path,wpm_eff) output_wpm_t50];
[tmp ind] = sort(sum(WPM_pv<=0.05,2),'descend');
output_wpm_t50 = output_wpm_t50(ind,:);
output_wpm_t50 = output_wpm_t50(1:nnz(tmp>=2),:);

writetable(output_bpm_t50,sprintf('%s/results_collection/yeast/yeast_across_pheno_pvalues_t50.xls',getenv('BRIDGEPATH')),'Sheet',1)
writetable(output_wpm_t50,sprintf('%s/results_collection/yeast/yeast_across_pheno_pvalues_t50.xls',getenv('BRIDGEPATH')),'Sheet',2)

clear BPM_tmp WPM_tmp PATH_tmp BPM_pv WPM_pv PATH_pv
