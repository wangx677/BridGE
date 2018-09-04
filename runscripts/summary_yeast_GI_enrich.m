phenos = {'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};
models = {'DD','RD','RR','combined'};

GI_enrichment_t25 = zeros(1,12);
GI_enrichment_t50 = zeros(1,12);

colnames = {'bpm_N_EE_neg_pv','bpm_N_EE_pos_pv','bpm_N_EN_neg_pv','bpm_N_EN_pos_pv','bpm_N_NN_neg_pv','bpm_N_NN_pos_pv','wpm_N_EE_neg_pv','wpm_N_EE_pos_pv','wpm_N_EN_neg_pv','wpm_N_EN_pos_pv','wpm_N_NN_neg_pv','wpm_N_NN_pos_pv'};

k=1;
for i=1:length(phenos)
     for j=1:length(models)
          files{k} = sprintf('%s_%s',phenos{i},models{j});
          projectdir=sprintf('%s/project_yeast_%s_complex_t25_b50_mhygeSSI',getenv('BRIDGEPATH'),phenos{i});
          data = readtable(sprintf('%s/output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',projectdir,models{j}),'Sheet',2);
          if isempty(data)~=1
               colnames_new = data.Row(find(ismember(data.Row,colnames)));
               ind1 = find(ismember(colnames,colnames_new));
               ind2 = find(ismember(data.Row,colnames_new));
               GI_enrichment_t25(k,ind1) = round(-log10(data.fdr25(ind2)')*100)/100;
          else
               GI_enrichment_t25(k,:) = 0;
          end
          

          projectdir=sprintf('%s/project_yeast_%s_complex_t50_b25_mhygeSSI',getenv('BRIDGEPATH'),phenos{i});
          data = readtable(sprintf('%s/output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',projectdir,models{j}),'Sheet',2);
          if isempty(data)~=1
               colnames_new = data.Row(find(ismember(data.Row,colnames)));
               ind1 = find(ismember(colnames,colnames_new));
               ind2 = find(ismember(data.Row,colnames_new));
               GI_enrichment_t50(k,ind1) = round(-log10(data.fdr25(ind2)')*100)/100;
          else
               GI_enrichment_t50(k,:) = 0;
          end
          k = k+1;
     end
end
     
GI_enrichment_t50(isinf(GI_enrichment_t50)) = 16;
GI_enrichment_t25(isinf(GI_enrichment_t25)) = 16;

files = files';

GI_enrichment_t25 = array2table(GI_enrichment_t25);
GI_enrichment_t25.Properties.VariableNames = colnames;
GI_enrichment_t25.Properties.RowNames = files;

GI_enrichment_t50 = array2table(GI_enrichment_t50);
GI_enrichment_t50.Properties.VariableNames = colnames;
GI_enrichment_t50.Properties.RowNames = files;

outputfile = sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_t25.xls',getenv('BRIDGEPATH'));
writetable(GI_enrichment_t25,outputfile,'WriteRowNames',true);

outputfile = sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_t50.xls',getenv('BRIDGEPATH'));
writetable(GI_enrichment_t50,outputfile,'WriteRowNames',true);

t25 = sum(table2array(GI_enrichment_t25)>=-log10(0.05));
t50 = sum(table2array(GI_enrichment_t50)>=-log10(0.05));

bpm_t25 = t25(1:6);wpm_t25 = t25(7:end);
bpm_t50 = t50(1:6);wpm_t50 = t50(7:end);

close all
bar([bpm_t25;bpm_t50]')
setfig
set(gca,'XtickLabels',{'EE-','EE+','EN-','EN+','NN-','NN+'});
xlabel('GI networks')
ylabel(sprintf('# of enriched results (p<0.05) out of %s',num2str(size(GI_enrichment_t25,1))));
legend({'t25-b50','t50-b25'})
saveas(gcf,sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_BPM.png',getenv('BRIDGEPATH')))
saveas(gcf,sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_BPM.pdf',getenv('BRIDGEPATH')))

close all
bar([wpm_t25;wpm_t50]')
setfig
set(gca,'XtickLabels',{'EE-','EE+','EN-','EN+','NN-','NN+'});
xlabel('GI networks')
ylabel(sprintf('# of enriched results (p<0.05) out of %s',num2str(size(GI_enrichment_t25,1))));
legend({'t25-b50','t50-b25'})
saveas(gcf,sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_WPM.png',getenv('BRIDGEPATH')))
saveas(gcf,sprintf('%s/results_collection/yeast/summary_yeast_GI_enrich_WPM.pdf',getenv('BRIDGEPATH')))

