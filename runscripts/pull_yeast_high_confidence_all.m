phenos = {'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};
models = {'DD','RD','RR','combined'};
binarization = {'t25_b50','t50_b25'};

% same pathway standard was used across data
for m = 1:length(binarization)
for i=1:length(phenos)
    projectdir=sprintf('%s/project_yeast_%s_complex_%s_mhygeSSI',getenv('BRIDGEPATH'),phenos{i},binarization{m});

    for j=1:length(models)
          try
               data_bpm = readtable(sprintf('%s/output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',projectdir,models{j}),'Sheet',3);
          catch
               data_bpm = '';
          end

          if isempty(data_bpm)~=1
               tmp = table2array(data_bpm(:,{'GI_enrich_pv_EE_neg','GI_enrich_pv_EE_pos','GI_enrich_pv_EN_neg','GI_enrich_pv_EN_pos','GI_enrich_pv_NN_neg','GI_enrich_pv_NN_pos'}));
               ind_bpm = find(min(tmp,[],2)<=0.05); % keep BPMs that can be enriched in SGA
               data_bpm = data_bpm(ind_bpm,:);
               enrichment = min(tmp(ind_bpm,:),[],2);
               enrichment = table(enrichment);
               diseasemodel = repmat(models(j),length(ind_bpm),1);
               diseasemodel = table(diseasemodel);
               data_bpm = data_bpm(:,{'path1','path2','bpm_size','fdrBPM','eff_bpm','bpm_pv_discovery','bpm_ranksum_discovery'});
               data_bpm.Properties.VariableNames = {'path1','path2','size','fdr','effect','p','ranksum'};
               data_bpm = [data_bpm enrichment diseasemodel];

          end

          try
               data_wpm = readtable(sprintf('%s/output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',projectdir,models{j}),'Sheet',4);
          catch
               data_wpm = '';
          end
     
          if isempty(data_wpm)~=1
               tmp = table2array(data_wpm(:,{'GI_enrich_pv_EE_neg','GI_enrich_pv_EE_pos','GI_enrich_pv_EN_neg','GI_enrich_pv_EN_pos','GI_enrich_pv_NN_neg','GI_enrich_pv_NN_pos'}));
               ind_wpm = find(min(tmp,[],2)<=0.05); % keep WPMs that can be enriched in SGA
               data_wpm = data_wpm(ind_wpm,:);
               enrichment = min(tmp(ind_wpm,:),[],2);
               enrichment = table(enrichment);
               diseasemodel = repmat(models(j),length(ind_wpm),1);
               diseasemodel = table(diseasemodel);
               data_wpm = data_wpm(:,{'path','wpm_size','fdrWPM','eff_wpm','wpm_pv_discovery','wpm_ranksum_discovery'});
               data_wpm.Properties.VariableNames = {'path1','size','fdr','effect','p','ranksum'};
               path2 = data_wpm.path1;
               path2 = table(path2);
               data_wpm = [data_wpm(:,1) path2 data_wpm(:,2:end) enrichment diseasemodel];
          end
     
          if isempty(data_bpm)~=1 & isempty(data_wpm)~=1
               data{j} = [data_bpm;data_wpm];
          elseif isempty(data_bpm)~=1
               data{j} = data_bpm;
          elseif isempty(data_wpm)~=1
               data{j} = data_wpm;
          else
               data{j} = '';
          end
          clear data_bpm data_wpm
     end

     data(find(cellfun(@(x)isempty(x),data)==1))=[];
     results_tmp = data{1};
     for k=2:length(data)
          results_tmp = [results_tmp;data{k}];
     end
    
     % add experimental condition 
     condition = repmat(phenos(i),size(results_tmp,1),1);
     condition = table(condition);
     results_tmp = [results_tmp condition];

     % keep FDR<0.25
     ind = find(results_tmp.fdr<=0.25);
     results_tmp = results_tmp(ind,:);

     % keep p<0.005
     ind = find(results_tmp.p<=0.005);
     results_tmp = results_tmp(ind,:);

     % sort BPMs by name and their stats
     results_tmp = sortrows(results_tmp,{'path1','path2','p','fdr'},{'ascend','ascend','ascend','ascend'});

     % keep unique discoveries
     tmp = cellfun(@(x,y)strjoin([{x},{y}],':'),results_tmp.path1,results_tmp.path2,'uniform',0);
     [tmp idx] = unique(tmp);
     results_tmp = results_tmp(idx,:);
     
     % sort BPMs by stats
     results_tmp = sortrows(results_tmp,{'p'},{'ascend'});
     results{i} = results_tmp;
%     writetable(results{i},sprintf('/project/csbio/wwang/BridGE/results_collection/yeast/results_high_confidence_%s.xls',binarization{m}),'Sheet',i);
     clear data results_tmp tmp idx
end

results_all_tmp = results{1};
for k=2:length(results)
     results_all_tmp = [results_all_tmp;results{k}];
end

results_all = sortrows(results_all_tmp,{'p'},{'ascend'});
writetable(results_all,sprintf('/project/csbio/wwang/BridGE/results_collection/yeast/results_high_confidence_%s.xls',binarization{m}));
writetable(results_all,sprintf('/project/csbio/wwang/BridGE/results_collection/yeast/results_high_confidence_%s.csv',binarization{m}),'Delimiter',';');
clear results results_all results_all_tmp
end
