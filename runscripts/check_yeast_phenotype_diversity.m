pheno={'SC4NQO01ugml_38h','SCCHX05ugml_38h','SCpH3_38h','SCpH8_38h','YPD42_40h','YPDCHX05_40h','YPDSDS_40h','YPGLYCEROL_40h'};
binary={'t25_b50','t50_b25','t25_b25'};
groups = readtable('/project/csbio/wwang/yeast_1011strains/maf005_groups.txt','FileType','text','ReadVariableNames',false);
groups.Var2(find(cellfun(@(x)isempty(x),table2array(groups(:,2)))))={'other'};
uniquegroups = unique(groups.Var2);

origin = readtable('/project/csbio/wwang/yeast_1011strains/1011_strains_info_forWen_withPloidyInfo.xlsx');
origin = origin(:,[6,9]);

for k=1:length(binary)
     %(1) check how frequently a yeast strain's growth ratio is in the extreme end across different conditions
     tb = split(binary{k},'_');
     t = split(tb{1},'t');
     t = t{2};
     b = split(tb{2},'b');
     b = b{2};

     for i = 1:length(pheno) 
          load(sprintf('project_yeast_%s_complex_%s_mhygeSSI/gwas_data_final.mat',pheno{i},binary{k}));
          strains{i} = SNPdata.fid;
          phenos{i} = SNPdata.pheno;
     end

     strains_all = union(strains{1},strains{2});
     for i=3:length(pheno)
          strains_all = union(strains_all,strains{i});
     end

     output = nan(length(strains_all),length(pheno));
     for i=1:length(pheno)
          ind = find(ismember(strains_all,strains{i}(find(phenos{i}==1))));
          output(ind,i) = 1;
          ind = find(ismember(strains_all,strains{i}(find(phenos{i}==0))));
          output(ind,i) = -1;
     end

     subplot(2,1,1)
     S = sum(output==1,2);
     x = unique(S);
     for i=1:length(x)
          y(i) = nnz(S==x(i));
     end
     bar(x,y)
     setfig
     ylabel('counts')
     xlabel('number of conditions')
     title(sprintf('top %s% growth ratio',t))
     clear x y

     subplot(2,1,2)
     S = sum(output==-1,2);
     x = unique(S);
     for i=1:length(x)
          y(i) = nnz(S==x(i));
     end
     bar(x,y)
     setfig
     ylabel('counts')
     xlabel('number of conditions')
     title(sprintf('bottom %s% growth ratio',b))
     tightfig
     clear x y

     saveas(gcf,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity1_%s.pdf',getenv('BRIDGEPATH'),binary{k}))
     saveas(gcf,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity1_%s.png',getenv('BRIDGEPATH'),binary{k}))
     close all

     % (2) put together of (1)
     freqmat = nan(8,8);
     for i=0:8
          for j=0:8
               freqmat(i+1,j+1) = nnz(sum(output==1,2)==i & sum(output==-1,2)==j);
          end
     end
     bar3(fliplr(freqmat))
     set(gca,'XTickLabel',[8:-1:0])
     set(gca,'YTickLabel',[0:8])
     setfig
     xlabel('number of conditions with high growth ratio)')
     ylabel('number of conditions with low growth ratio)')
     zlabel('counts')
     title('frequency of conditions yeast strains are in high/low end')
     tightfig
     saveas(gcf,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity2_%s.pdf',getenv('BRIDGEPATH'),binary{k}))
     saveas(gcf,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity2_%s.png',getenv('BRIDGEPATH'),binary{k}))
     close all

     % (3) pairwise comparison of conditions
     phenosim = zeros(size(output,2),size(output,2));
     for i=1:(size(output,2)-1);
          for j=i+1:size(output,2)
               phenosim(i,j) =  nnz(output(:,i)-output(:,j)==0)/max(sum(isnan(output)~=1));
          end
     end
     phenosim = phenosim+phenosim';

     phenosim = round(phenosim*100)/100;
     phenosim = array2table(phenosim);
     phenosim.Properties.VariableNames = pheno;
     phenosim.Properties.RowNames = pheno;
     writetable(phenosim,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity3_%s.xls',getenv('BRIDGEPATH'),binary{k}),'WriteRowNames',1)

     % (4) check distributions of yeast groups in case and control
     for i=1:length(pheno)
          for m=1:length(uniquegroups)
              tmp = intersect(groups.Var1(find(ismember(groups.Var2,uniquegroups{m}))),strains{i});;
              casefreq(m) = length(intersect(tmp,strains{i}(find(phenos{i}==0))))/length(tmp);
              controlfreq(m) = length(intersect(tmp,strains{i}(find(phenos{i}==1))))/length(tmp);
              samplesize(m) = length(tmp);
          end
          samplesize(find(samplesize==0)) = 0.01;
          casefreq(isnan(casefreq))=0;
          controlfreq(isnan(controlfreq))=0;
          scatter(casefreq,controlfreq,samplesize,1:length(uniquegroups));
          xlim([0 1]);
          ylim([0 1]);
          setfig
          xlabel('% in low group')
          ylabel('% in high group')
          title(pheno{i},'Interpreter', 'none')
          tightfig
          saveas(gcf,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity4_%s_%s.pdf',getenv('BRIDGEPATH'),binary{k},pheno{i}))
          saveas(gcf,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity4_%s_%s.png',getenv('BRIDGEPATH'),binary{k},pheno{i}))
          close all
     end
     clearvars -except pheno binary groups uniquegroups casefreq controlfreq samplesize
end

% (5) phenotypic diversity based on quantitative traits
data = readtable('/project/csbio/wwang/yeast_1011strains/forWen/pheno1002.txt','FileType','text');
for i=1:length(pheno)
     tmp = strsplit(pheno{i},'_');
     condition = tmp{1};
     ctime = tmp{2};
     strains{i} = data.standardized_name(find(ismember(data.condition,condition) & ismember(data.time,ctime)));
     phenos{i} = data.mean_ratio(find(ismember(data.condition,condition) & ismember(data.time,ctime)));;
     test = cell2mat(cellfun(@(x)strcmp(x,'NA'),phenos{i},'uniform',0));
     if nnz(test)>0
          strains{i}(find(test==1)) = [];
          phenos{i}(find(test==1)) = [];
     end
end

strains_all = union(strains{1},strains{2});
for i=3:length(pheno)
     strains_all = union(strains_all,strains{i});
end

output = nan(length(strains_all),length(pheno));
for i=1:length(pheno)
     [tmp ind] = ismember(strains{i},strains_all);
     output(ind,i) = cell2mat(cellfun(@(x)str2num(x),phenos{i},'uniform',0));
end

strains_all = strains_all(find(sum(isnan(output),2)==0));
output =output(find(sum(isnan(output),2)==0),:);

% Spearman Rho (rank correlation)
corr1 = corr(output,'type','Spearman');
corr1 = tril(corr1,-1) + triu(corr1,1);
corr1 = round(corr1*100)/100;
corr1 = array2table(corr1);
corr1.Properties.VariableNames = pheno;
corr1.Properties.RowNames = pheno;
writetable(corr1,sprintf('%s/results_collection/yeast/check_yeast_phenotype_diversity5_allstrains.xls',getenv('BRIDGEPATH')),'WriteRowNames',1)

[tmp ind] = ismember(strains_all,groups.Var1);
groups_all = groups.Var2(ind);
for i=1:length(uniquegroups)
     ind = find(ismember(groups_all,uniquegroups{i}));
     groupsize(i) = length(ind);
     corr2{i} = corr(output(ind,:)');
end

corr2 = corr(zscore(output)');
cg = clustergram(corr2);
plot(cg)
setfig
saveas(gcf,sprintf('%s/results_collection/yeast/yeast_phenotype_clustergram.pdf',getenv('BRIDGEPATH')))
saveas(gcf,sprintf('%s/results_collection/yeast/yeast_phenotype_clustergram.png',getenv('BRIDGEPATH')))
rowlabels = cellfun(@(x)str2num(x),cg.RowLabels);
     
