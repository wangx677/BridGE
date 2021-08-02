dirs = {'project_PD_NGRC_pbody_0kb','project_PD_NGRC_pbody_50kb','project_PD_NGRC_pbody_merge_0kb','project_PD_NGRC_pbody_merge_50kb','project_PD_Simon_pbody_0kb','project_PD_Simon_pbody_50kb','project_PD_Simon_pbody_merge_0kb','project_PD_Simon_pbody_merge_50kb'};

names = {'NGRC_0kb','NGRC_50kb','NGRC_merge_0kb','NGRC_merge_50kb','NIA2_0kb','NIA2_50kb','NIA2_merge_0kb','NIA2_merge_50kb'};

model = {'RR','DD','RD','combined'};

% summarize the results based on FDR<=0.25
for i=1:length(dirs)
     load(sprintf('/project/csbio/wwang/BridGE/%s/BPMind.mat',dirs{i}));
     for j=1:length(model)
          load(sprintf('/project/csbio/wwang/BridGE/%s/results_pbody_%s.mat',dirs{i},model{j}));
          for tt=1:2
               fdr_summary{tt}(j,i) = nnz(fdr_combined{tt}<0.255);
          end
     end
end

fdr_summary_protective = array2table(fdr_summary{1},'RowNames',model,'VariableNames',names);
fdr_summary_risk = array2table(fdr_summary{2},'RowNames',model,'VariableNames',names);           
cd('/project/csbio/wwang/BridGE/scripts_PD_pbody');
writetable(fdr_summary_protective,'pbody_fdr_summary.xls','WriteRowNames',true,'Sheet',1)
writetable(fdr_summary_risk,'pbody_fdr_summary.xls','WriteRowNames',true,'Sheet',2)

% pull all interactions with FDR<0.25 (risk only)
load /project/csbio/wwang/BridGE/project_PD_NGRC_pbody_50kb/BPMind.mat
genes1 = [WPM.pathway(BPM.path1idx);WPM.pathway];
genes2 = [WPM.pathway(BPM.path2idx);WPM.pathway];
clear BPM WPM

genes1_genes2 = cellfun(@(x,y)sprintf('%s_%s',x,y),genes1,genes2,'uniform',0);

pbody_GI_table = nan(length(genes1),24);
model_table = repmat({''},length(genes1),length(dirs));

for i=1:length(dirs)
     load(sprintf('/project/csbio/wwang/BridGE/%s/BPMind.mat',dirs{i}));
     tmp_genes1 = [WPM.pathway(BPM.path1idx);WPM.pathway];
     tmp_genes2 = [WPM.pathway(BPM.path2idx);WPM.pathway]; 
     tmp_genes1_genes2 = cellfun(@(x,y)sprintf('%s_%s',x,y),tmp_genes1,tmp_genes2,'uniform',0);
    
     for j=1:length(model)
          load(sprintf('/project/csbio/wwang/BridGE/%s/results_pbody_%s.mat',dirs{i},model{j}));
          tmp_density(:,j) = [BPM_density{2} WPM_density{2}]';
          tmp_pv(:,j) = [BPM_density_pv{2} WPM_density_pv{2}]';    
          tmp_fdr(:,j) = fdr_combined{2}';
     end
     
     [tmpfdr,idx]=min(tmp_fdr,[],2);
     
     for k=1:length(idx)
          maxdensity(k) = tmp_density(k,idx(k));
          minpv(k) = tmp_pv(k,idx(k));
          minfdr(k) = tmp_fdr(k,idx(k));
          minmodel{k} = model{idx(k)};
     end

     [tmp idxa idxb] = intersect(genes1_genes2,tmp_genes1_genes2);
     pbody_GI_table(idxa,3*(i-1)+1) = round(maxdensity(idxb),2);
     pbody_GI_table(idxa,3*(i-1)+2) = minpv(idxb);
     pbody_GI_table(idxa,3*(i-1)+3) = round(minfdr(idxb),2);

     model_table(idxa,i) = minmodel(idxb);
     clear maxdensity minpv minfdr minmodel tmp*
end

pbody_GI_table = [array2table(genes1) array2table(genes2) array2table(pbody_GI_table,'VariableNames',{'NGRC_0kb_density','NGRC_0kb_pv','NGRC_0kb_FDR','NGRC_50kb_density','NGRC_50kb_pv','NGRC_50kb_FDR','NGRC_merge_0kb_density','NGRC_merge_0kb_pv','NGRC_merge_0kb_FDR','NGRC_merge_50kb_density','NGRC_merge_50kb_pv','NGRC_merge_50kb_FDR','NIA2_0kb_density','NIA2_0kb_pv','NIA2_0kb_FDR','NIA2_50kb_density','NIA2_50kb_pv','NIA2_50kb_FDR','NIA2_merge_0kb_density','NIA2_merge_0kb_pv','NIA2_merge_0kb_FDR','NIA2_merge_50kb_density','NIA2_merge_50kb_pv','NIA2_merge_50kb_FDR'})];


% filtering
NGRC_pv_sum = sum([pbody_GI_table.NGRC_0kb_pv pbody_GI_table.NGRC_50kb_pv pbody_GI_table.NGRC_merge_0kb_pv pbody_GI_table.NGRC_merge_50kb_pv]<=0.05,2);
NIA2_pv_sum = sum([pbody_GI_table.NIA2_0kb_pv pbody_GI_table.NIA2_50kb_pv pbody_GI_table.NIA2_merge_0kb_pv pbody_GI_table.NIA2_merge_50kb_pv]<=0.05,2);

NGRC_FDR_sum = sum([pbody_GI_table.NGRC_0kb_FDR pbody_GI_table.NGRC_50kb_FDR pbody_GI_table.NGRC_merge_0kb_FDR pbody_GI_table.NGRC_merge_50kb_FDR]<0.255,2);
NIA2_FDR_sum = sum([pbody_GI_table.NIA2_0kb_FDR pbody_GI_table.NIA2_50kb_FDR pbody_GI_table.NIA2_merge_0kb_FDR pbody_GI_table.NIA2_merge_50kb_FDR]<0.255,2);

validation = ((NGRC_FDR_sum>=1 & NIA2_pv_sum>=1) | (NIA2_FDR_sum>=1 & NGRC_pv_sum>=1));

ind = find(NGRC_pv_sum>=1 | NIA2_pv_sum>=1);

pbody_GI_table = [pbody_GI_table table(NGRC_pv_sum) table(NIA2_pv_sum) table(NGRC_FDR_sum) table(NIA2_FDR_sum) table(validation)];

pbody_GI_table_filter_pv = pbody_GI_table(ind,:);
model_table = array2table(model_table,'VariableNames',names);
model_table_filter_pv = model_table(ind,:);

ind = find(NGRC_FDR_sum>=1 | NIA2_FDR_sum>=1);

pbody_GI_table_filter_FDR = pbody_GI_table(ind,:);
model_table_filter_FDR = model_table(ind,:);

save generate_results pbody_GI_table pbody_GI_table_filter* model_table model_table_filter*

writetable(pbody_GI_table_filter_pv,'pbody_GI_table_filter_pv.xls')
writetable(pbody_GI_table_filter_FDR,'pbody_GI_table_filter_FDR.xls')
