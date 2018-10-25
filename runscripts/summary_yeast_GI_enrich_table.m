direction = 'same'; % this needs to be changed to something else if doesn't require to match interaction direction

phenos = {'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};
models = {'DD','RD','RR','combined'};

GI_enrichment = nan(12,24); % match the format in slide "Summary of complex interactions"
GI_enrichment = array2table(GI_enrichment);
GI_enrichment.Properties.RowNames={'DD_bpm','DD_wpm','DD_path','RD_bpm','RD_wpm','RD_path','RR_bmp','RR_wpm','RR_path','combined_bpm','combined_wpm','combined_paht'};

k = 1;
for i = 1:length(phenos)
     for j = {'t25b50','t25b25','t50b25'}
          GI_enrichment.Properties.VariableNames{k} = sprintf('%s_%s',phenos{i},j{1});
          k = k+1;
     end
end

bpmcolnames = {'bpm_N_EE_neg_pv','bpm_N_EE_pos_pv','bpm_N_EN_neg_pv','bpm_N_EN_pos_pv','bpm_N_NN_neg_pv','bpm_N_NN_pos_pv'};
wpmcolnames = {'wpm_N_EE_neg_pv','wpm_N_EE_pos_pv','wpm_N_EN_neg_pv','wpm_N_EN_pos_pv','wpm_N_NN_neg_pv','wpm_N_NN_pos_pv'};

k=1;
for i=1:length(phenos)
     for j=1:length(models)
          files{k} = sprintf('%s_%s',phenos{i},models{j});
          projectdir=sprintf('%s/project_yeast_%s_complex_t25_b50_mhygeSSI',getenv('BRIDGEPATH'),phenos{i});
          if strcmp(direction,'same')~=1
               data1 = readtable(sprintf('%s/output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',projectdir,models{j}),'Sheet',2);
          else
               filename = sprintf('%s/GI_enrichment_samedirection_%s.xls',projectdir,models{j});
               if exist(filename,'file')==2
                    data1 = readtable(sprintf('%s/GI_enrichment_samedirection_%s.xls',projectdir,models{j}));
               else
                    data1 = '';
               end
          end

          projectdir=sprintf('%s/project_yeast_%s_complex_t25_b25_mhygeSSI',getenv('BRIDGEPATH'),phenos{i});
          if strcmp(direction,'same')~=1
               data2 = readtable(sprintf('%s/output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',projectdir,models{j}),'Sheet',2);
          else
               filename = sprintf('%s/GI_enrichment_samedirection_%s.xls',projectdir,models{j});
               if exist(filename,'file')==2
                    data2 = readtable(sprintf('%s/GI_enrichment_samedirection_%s.xls',projectdir,models{j}));
               else
                    data2 = '';
               end

          end

          projectdir=sprintf('%s/project_yeast_%s_complex_t50_b25_mhygeSSI',getenv('BRIDGEPATH'),phenos{i});
          if strcmp(direction,'same')~=1
               data3 = readtable(sprintf('%s/output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',projectdir,models{j}),'Sheet',2);
          else
               filename = sprintf('%s/GI_enrichment_samedirection_%s.xls',projectdir,models{j});
               if exist(filename,'file')==2
                    data3 = readtable(sprintf('%s/GI_enrichment_samedirection_%s.xls',projectdir,models{j}));
               else
                    data3 = '';
               end
          end

          if isempty(data1)~=1
               bpmcolnames_new = data1.Row(find(ismember(data1.Row,bpmcolnames)));
               ind1 = find(ismember(bpmcolnames,bpmcolnames_new));
               ind2 = find(ismember(data1.Row,bpmcolnames_new));
               if isempty(ind1)~=1
                    GI_enrichment((j-1)*3+1,(i-1)*3+1) = {max(round(-log10(data1.fdr25(ind2)')*100)/100)>=-log10(0.05)};
               else
                    GI_enrichment((j-1)*3+1,(i-1)*3+1) = {0};
               end

               wpmcolnames_new = data1.Row(find(ismember(data1.Row,wpmcolnames)));
               ind1 = find(ismember(wpmcolnames,wpmcolnames_new));
               ind2 = find(ismember(data1.Row,wpmcolnames_new));
               if isempty(ind1)~=1
                    GI_enrichment((j-1)*3+2,(i-1)*3+1) = {max(round(-log10(data1.fdr25(ind2)')*100)/100)>=-log10(0.05)};
               else
                    GI_enrichment((j-1)*3+2,(i-1)*3+1) = {0};
               end
          else
               GI_enrichment((j-1)*3+1,(i-1)*3+1) = {0};
               GI_enrichment((j-1)*3+2,(i-1)*3+1) = {0};
          end
     
          if isempty(data2)~=1
               bpmcolnames_new = data2.Row(find(ismember(data2.Row,bpmcolnames)));
               ind1 = find(ismember(bpmcolnames,bpmcolnames_new));
               ind2 = find(ismember(data2.Row,bpmcolnames_new));
               if isempty(ind1)~=1
                    GI_enrichment((j-1)*3+1,(i-1)*3+2) = {max(round(-log10(data2.fdr25(ind2)')*100)/100)>=-log10(0.05)};
               else
                    GI_enrichment((j-1)*3+1,(i-1)*3+2) = {0};
               end

               wpmcolnames_new = data2.Row(find(ismember(data2.Row,wpmcolnames)));
               ind1 = find(ismember(wpmcolnames,wpmcolnames_new));
               ind2 = find(ismember(data2.Row,wpmcolnames_new));
               if isempty(ind1)~=1
                    GI_enrichment((j-1)*3+2,(i-1)*3+2) = {max(round(-log10(data2.fdr25(ind2)')*100)/100)>=-log10(0.05)};
               else
                    GI_enrichment((j-1)*3+2,(i-1)*3+2) = {0};
               end
          else
               GI_enrichment((j-1)*3+1,(i-1)*3+2) = {0};
               GI_enrichment((j-1)*3+2,(i-1)*3+2) = {0};
          end

          if isempty(data3)~=1
               bpmcolnames_new = data3.Row(find(ismember(data3.Row,bpmcolnames)));
               ind1 = find(ismember(bpmcolnames,bpmcolnames_new));
               ind2 = find(ismember(data3.Row,bpmcolnames_new));
               if isempty(ind1)~=1
                    GI_enrichment((j-1)*3+1,(i-1)*3+3) = {max(round(-log10(data3.fdr25(ind2)')*100)/100)>=-log10(0.05)};
               else
                    GI_enrichment((j-1)*3+1,(i-1)*3+3) = {0};
               end

               wpmcolnames_new = data3.Row(find(ismember(data3.Row,wpmcolnames)));
               ind1 = find(ismember(wpmcolnames,wpmcolnames_new));
               ind2 = find(ismember(data3.Row,wpmcolnames_new));
               if isempty(ind1)~=1
                    GI_enrichment((j-1)*3+2,(i-1)*3+3) = {max(round(-log10(data3.fdr25(ind2)')*100)/100)>=-log10(0.05)};
               else
                    GI_enrichment((j-1)*3+2,(i-1)*3+3) = {0};
               end
          else
               GI_enrichment((j-1)*3+1,(i-1)*3+3) = {0};
               GI_enrichment((j-1)*3+2,(i-1)*3+3) = {0};
          end
     end
     k = k+1;
end
     
if strcmp(direction,'same')~=1
     outputfile = sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_table.xls',getenv('BRIDGEPATH'));
     writetable(GI_enrichment,outputfile,'WriteRowNames',true);
else
     outputfile = sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_table_same_direction.xls',getenv('BRIDGEPATH'));
     writetable(GI_enrichment,outputfile,'WriteRowNames',true);
end
