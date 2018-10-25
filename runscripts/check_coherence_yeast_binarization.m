binary1 = 't50_b25';
binary2 = 't25_b25';

phenos = {'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};
models = {'DD','RD','RR','combined'};

j=1;
for i=1:length(phenos)
     for k=1:length(models)
          discoverydir = sprintf('/project/csbio/wwang/BridGE/project_yeast_%s_complex_%s_mhygeSSI',phenos{i},binary1);
          validationdir = sprintf('/project/csbio/wwang/BridGE/project_yeast_%s_complex_%s_mhygeSSI',phenos{i},binary2);
          resultfile = sprintf('results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',models{k});
          [logp1, No_total1, No_discovery1, No_valid1, No_overlap1] = test_agreement_across_pathway_interaction(discoverydir,validationdir,resultfile,'snp_pathway_min5_max300.mat','BPMind.mat','top');
          [logp2, No_total2, No_discovery2, No_valid2, No_overlap2] = test_agreement_across_pathway_interaction(discoverydir,validationdir,resultfile,'snp_pathway_min5_max300.mat','BPMind.mat','fdr');
          
          files{j} = sprintf('%s_%s_BPM',phenos{i},models{k}); 
          output_top{j} = sprintf('%.2f(t25:%d; t50:%d; overlap:%d)',logp1.BPM,No_discovery1.BPM,No_valid1.BPM,No_overlap1.BPM);
          output_fdr{j} = sprintf('%.2f(t25:%d; t50:%d; overlap:%d)',logp2.BPM,No_discovery2.BPM,No_valid2.BPM,No_overlap2.BPM);
          j = j+1;

          files{j} = sprintf('%s_%s_WPM',phenos{i},models{k});
          output_top{j} = sprintf('%.2f(t25:%d; t50:%d; overlap:%d)',logp1.WPM,No_discovery1.WPM,No_valid1.WPM,No_overlap1.WPM);
          output_fdr{j} = sprintf('%.2f(t25:%d; t50:%d; overlap:%d)',logp2.WPM,No_discovery2.WPM,No_valid2.WPM,No_overlap2.WPM);
          j = j+1;

          files{j} = sprintf('%s_%s_PATH',phenos{i},models{k});
          output_top{j} = sprintf('%.2f(t25:%d; t50:%d; overlap:%d)',logp1.PATH,No_discovery1.PATH,No_valid1.PATH,No_overlap1.PATH);
          output_fdr{j} = sprintf('%.2f(t25:%d; t50:%d; overlap:%d)',logp2.PATH,No_discovery2.PATH,No_valid2.PATH,No_overlap2.PATH);
          j = j+1;
     end
end

files = files';
output_top = output_top';
output_fdr = output_fdr';
output = table(files,output_top,output_fdr);
writetable(output,sprintf('%s/results_collection/yeast/yeast_coherence_%s_vs_%s.xls',getenv('BRIDGEPATH'),binary1,binary2))


output_top_logp = zeros(1,length(output_top));
output_top_overlap = zeros(1,length(output_top));
for i=1:length(output_top)
     tmp = strsplit(output_top{i},'(');
     output_top_logp(i) = str2num(tmp{1});
     tmp = strsplit(tmp{2},':');
     tmp = strsplit(tmp{end},')');
     output_top_overlap(i) = str2num(tmp{1});
     tmp = strsplit(files{i},'_');
     interaction{i} = tmp{end};
end

close all
subplot(2,2,1)
hist(output_top_overlap(find(ismember(interaction,'BPM'))))
setfig
xlabel('overlap distribution (out of 100)')
ylabel('frequency')
title('BPM')

subplot(2,2,2)
hist(output_top_overlap(find(ismember(interaction,'WPM'))))
setfig
xlabel('overlap distribution (out of 20)')
ylabel('frequency')
title('WPM')

subplot(2,2,3)
hist(output_top_logp(find(ismember(interaction,'BPM'))))
setfig
xlabel('-log10 pv distribution (hypergeometric test)')
ylabel('frequency')
title('BPM')

subplot(2,2,4)
hist(output_top_logp(find(ismember(interaction,'WPM'))))
setfig
xlabel('-log10 pv distribution (hypergeometric test)')
ylabel('frequency')
title('WPM')
saveas(gcf,sprintf('%s/results_collection/yeast/yeast_coherence_%s_vs_%s.png',getenv('BRIDGEPATH'),binary1,binary2))
saveas(gcf,sprintf('%s/results_collection/yeast/yeast_coherence_%s_vs_%s.pdf',getenv('BRIDGEPATH'),binary1,binary2))


