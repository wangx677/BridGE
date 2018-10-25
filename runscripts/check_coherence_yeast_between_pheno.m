phenos = {'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};
models = {'DD','RD','RR','combined'};
binary = {'t25_b50','t50_b25','t25_b25'};

for n=1:length(binary)
     for i=1:length(phenos)-1
          for m=i+1:length(phenos)
               for k=1:length(models)
                    discoverydir = sprintf('/project/csbio/wwang/BridGE/project_yeast_%s_complex_%s_mhygeSSI',phenos{i},binary{n});
                    validationdir = sprintf('/project/csbio/wwang/BridGE/project_yeast_%s_complex_%s_mhygeSSI',phenos{m},binary{n});
                    resultfile = sprintf('results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',models{k});
                    [logp1, No_total1, No_discovery1, No_valid1, No_overlap1] = test_agreement_across_pathway_interaction(discoverydir,validationdir,resultfile,'snp_pathway_min5_max300.mat','BPMind.mat','top');
                    [logp2, No_total2, No_discovery2, No_valid2, No_overlap2] = test_agreement_across_pathway_interaction(discoverydir,validationdir,resultfile,'snp_pathway_min5_max300.mat','BPMind.mat','fdr');
          
                    output_top_BPM{k}{i,m} = sprintf('%.2f(%d/%d)',logp1.BPM,No_overlap1.BPM,No_discovery1.BPM);
                    output_fdr_BPM{k}{i,m} = sprintf('%.2f(p1:%d; p2:%d; o:%d)',logp2.BPM,No_discovery2.BPM,No_valid2.BPM,No_overlap2.BPM);

                    output_top_WPM{k}{i,m} = sprintf('%.2f(%d/%d)',logp1.WPM,No_overlap1.WPM,No_discovery1.WPM);
                    output_fdr_WPM{k}{i,m} = sprintf('%.2f(p1:%d; p2:%d; o:%d)',logp2.WPM,No_discovery2.WPM,No_valid2.WPM,No_overlap2.WPM);

                    output_top_PATH{k}{i,m} = sprintf('%.2f(%d/%d)',logp1.PATH,No_overlap1.PATH,No_discovery1.PATH);
                    output_fdr_PATH{k}{i,m} = sprintf('%.2f(p1:%d; p1:%d; o:%d)',logp2.PATH,No_discovery2.PATH,No_valid2.PATH,No_overlap2.PATH);
               end
          end
     end

     for k=1:length(models)
          output_top_BPM{k}(length(phenos),:) = output_top_BPM{k}(1,1);
          output_fdr_BPM{k}(length(phenos),:) = output_fdr_BPM{k}(1,1);
          output_top_WPM{k}(length(phenos),:) = output_top_WPM{k}(1,1);
          output_fdr_WPM{k}(length(phenos),:) = output_fdr_WPM{k}(1,1);
          output_top_PATH{k}(length(phenos),:) = output_top_PATH{k}(1,1);
          output_fdr_PATH{k}(length(phenos),:) = output_fdr_PATH{k}(1,1);
     end

     output1 = array2table([output_top_BPM{1}; output_top_WPM{1}; output_top_PATH{1}]);
     output1.Properties.VariableNames=phenos;
     output1.Properties.RowNames = [cellfun(@(x)sprintf('BPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('WPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('PATH_%s',x),phenos,'uniform',0)];

     writetable(output1,sprintf('%s/results_collection/yeast/yeast_coherence_different_pheno_%s.xls',getenv('BRIDGEPATH'),binary{n}),'Sheet',1,'WriteRowNames',true);

     output1 = array2table([output_top_BPM{2}; output_top_WPM{2}; output_top_PATH{2}]);
     output1.Properties.VariableNames=phenos;
     output1.Properties.RowNames = [cellfun(@(x)sprintf('BPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('WPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('PATH_%s',x),phenos,'uniform',0)];

     writetable(output1,sprintf('%s/results_collection/yeast/yeast_coherence_different_pheno_%s.xls',getenv('BRIDGEPATH'),binary{n}),'Sheet',2,'WriteRowNames',true);

     output1 = array2table([output_top_BPM{3}; output_top_WPM{3}; output_top_PATH{3}]);
     output1.Properties.VariableNames=phenos;
     output1.Properties.RowNames = [cellfun(@(x)sprintf('BPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('WPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('PATH_%s',x),phenos,'uniform',0)];

     writetable(output1,sprintf('%s/results_collection/yeast/yeast_coherence_different_pheno_%s.xls',getenv('BRIDGEPATH'),binary{n}),'Sheet',3,'WriteRowNames',true);

     output1 = array2table([output_top_BPM{4}; output_top_WPM{4}; output_top_PATH{4}]);
     output1.Properties.VariableNames=phenos;
     output1.Properties.RowNames = [cellfun(@(x)sprintf('BPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('WPM_%s',x),phenos,'uniform',0) cellfun(@(x)sprintf('PATH_%s',x),phenos,'uniform',0)];

     writetable(output1,sprintf('%s/results_collection/yeast/yeast_coherence_different_pheno_%s.xls',getenv('BRIDGEPATH'),binary{n}),'Sheet',4,'WriteRowNames',true);

     save(sprintf('%s/results_collection/yeast/yeast_coherence_different_pheno_%s.mat',getenv('BRIDGEPATH'),binary{n}),'output*')

     % plot result 

     for k=1:length(output_top_BPM)
          output_top_BPM_logp{k} = zeros(size(output_top_BPM{k},1),size(output_top_BPM{k},2));
          output_top_BPM_overlap{k} = output_top_BPM_logp{k};
          for i=1:size(output_top_BPM{k},1)-1
               for j=i+1:size(output_top_BPM{k},2)
                    tmp = strsplit(output_top_BPM{k}{i,j},'(');
                    output_top_BPM_logp{k}(i,j) = str2num(tmp{1});
                    tmp = strsplit(tmp{2},'/');
                    output_top_BPM_overlap{k}(i,j) = str2num(tmp{1});
               end
          end
          output_top_BPM_logp{k} = squareform(output_top_BPM_logp{k}');
          output_top_BPM_overlap{k} = squareform(output_top_BPM_overlap{k}');
     end

     for k=1:length(output_top_WPM)
          output_top_WPM_logp{k} = zeros(size(output_top_WPM{k},1),size(output_top_WPM{k},2));
          output_top_WPM_overlap{k} = output_top_WPM_logp{k};
          for i=1:size(output_top_WPM{k},1)-1
               for j=i+1:size(output_top_WPM{k},2)
                    tmp = strsplit(output_top_WPM{k}{i,j},'(');
                    output_top_WPM_logp{k}(i,j) = str2num(tmp{1});
                    tmp = strsplit(tmp{2},'/');
                    output_top_WPM_overlap{k}(i,j) = str2num(tmp{1});
               end
          end
          output_top_WPM_logp{k} = squareform(output_top_WPM_logp{k}');
          output_top_WPM_overlap{k} = squareform(output_top_WPM_overlap{k}');
     end

     output_top_BPM_logp = cell2mat(output_top_BPM_logp);
     output_top_BPM_overlap = cell2mat(output_top_BPM_overlap);

     output_top_WPM_logp = cell2mat(output_top_WPM_logp);
     output_top_WPM_overlap = cell2mat(output_top_WPM_overlap);

     close all
     subplot(2,2,1)
     hist(output_top_BPM_overlap)
     setfig
     xlabel('overlap distribution (out of 100)')
     ylabel('frequency')
     title('BPM')

     subplot(2,2,2)
     hist(output_top_WPM_overlap)
     setfig
     xlabel('overlap distribution (out of 20)')
     ylabel('frequency')
     title('WPM')

     subplot(2,2,3)
     hist(output_top_BPM_logp)
     setfig
     xlabel('-log10 pv distribution (hypergeometric test)')
     ylabel('frequency')
     title('BPM')

     subplot(2,2,4)
     hist(output_top_WPM_logp)
     setfig
     xlabel('-log10 pv distribution (hypergeometric test)')
     ylabel('frequency')
     title('WPM')

     tightfig
     saveas(gcf,sprintf('%s/results_collection/yeast/yeast_coherence_different_pheno_%s.png',getenv('BRIDGEPATH'),binary{n}))
     saveas(gcf,sprintf('%s/results_collection/yeast/yeast_coherence_different_pheno_%s.pdf',getenv('BRIDGEPATH'),binary{n}))
     close all
     clearvars -except phenos binary models n
end
