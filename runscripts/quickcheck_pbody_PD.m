cd /project/csbio/wwang/BridGE/project_PD_Simon_mhygeSSI

pbody={'LSM14A','EDC4','EDC3','DCP1','DCP2','PATIL','PNRC2','XRN1','PAN2','PAN3','LSM1','LSM2','LSM3','LSM4','LSM5','LSM6','LSM7','4E-T'};

load snpgenemapping_50kb.mat

%% get snp index that can be mapped to pbody
ind = find(ismember(snp2gene.genelist,pbody));
ind1 = find(sum(snp2gene.sgmatrix(:,ind),2)~=0);

%% find existing pathways that are signficantly overlapped with pbody SNPs
load snp_pathway_min10_max300.mat
a = size(snpset.spmatrix,1);
b = sum(snpset.spmatrix);
c = sum(snpset.spmatrix(ind1,:));
d = length(ind1);
pbody_enrich = hygetest(a,b,c,d);
ind2 = find(pbody_enrich>=-log10(0.05/length(snpset.pathwaynames)) & c>=10);
pbody_pathways = snpset.pathwaynames(ind2);

%% a table about pbody realted pathway
pbody_enrich = pbody_enrich(ind2);
pbody_snps = c(ind2);
pathway_snps = b(ind2);
for i=1:length(pbody_pathways)
     pbody_genes{i} = strjoin(intersect(snp2gene.genelist(find(sum(snp2gene.sgmatrix(find(snpset.spmatrix(:,find(ismember(snpset.pathwaynames,pbody_pathways{i})))==1),:))>0)),pbody),';');
end
pbody_genes = pbody_genes';
pbody_snps = pbody_snps';
pathway_snps = pathway_snps';
pbody_enrich = pbody_enrich';
pbody_pathway_table = table(pbody_pathways,pbody_genes,pbody_snps,pathway_snps,pbody_enrich);

load BPMind.mat
model = {'RR','DD','RD','combined'};
dirs = {'Simon','NGRC'};

%% check current BridGE results to any discoveries are related to pbody
% didn't find any discoveries (FDR<0.4) have pbody pathways

output_bpm = cell2table(cell(0,13),'VariableNames',{'data','disease_model','path1','path2','fdrBPM','eff_bpm','bpm_size','bpm_pv_discovery','bpm_ranksum_discovery','bpm_pv_valid','bpm_ranksum_valid','bpm_path1_drivers','bpm_path2_drivers'}); 
output_wpm = cell2table(cell(0,11),'VariableNames',{'data','disease_model','path','fdrWPM','eff_wpm','wpm_size','wpm_pv_discovery','wpm_ranksum_discovery','wpm_pv_valid','wpm_ranksum_valid','wpm_path_drivers'});
output_path = cell2table(cell(0,10),'VariableNames',{'data','disease_model','path','fdrPATH','eff_path','path_size','path_pv_discovery','path_ranksum_discovery','path_pv_valid','path_ranksum_valid'});

for j=1:length(dirs)
     cd(sprintf('/project/csbio/wwang/BridGE/project_PD_%s_mhygeSSI',dirs{j}))
     for i = 1:length(model)
          load(sprintf('output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',model{i}))
          if exist('output_bpm_table','var')
               ind_bpm = find(ismember(output_bpm_table.path1,pbody_pathways)|ismember(output_bpm_table.path2,pbody_pathways));
               if isempty(ind_bpm)==0
                    disease_model = repmat(model(i),length(ind_bpm),1);
                    data = repmat(dirs(j),length(ind_bpm),1);
                    output_bpm = [output_bpm;[table(data,disease_model),output_bpm_table(ind_bpm,:)]];
               end
          end

          if exist('output_wpm_table','var')
               ind_wpm = find(ismember(output_wpm_table.path,pbody_pathways));
               if isempty(ind_wpm)==0
                    disease_model = repmat(model(i),length(ind_wpm),1);
                    data = repmat(dirs(j),length(ind_wpm),1);
                    output_wpm = [output_wpm;[table(data,disease_model), output_wpm_table(ind_wpm,:)]];
               end
          end

          if exist('output_path_table','var')
               ind_path = find(ismember(output_path_table.path,pbody_pathways));
               if isempty(ind_path)==0
                    disease_model = repmat(model(i),length(ind_path),1);
                    data = repmat(dirs(j),length(ind_path),1);
                    output_path = [output_path; [table(data,disease_model), output_path_table(ind_path,:)]];
               end
          end

          clear output_bpm_table output_wpm_table output_path_table ind_bpm ind_wpm ind_path
     end
end

% relax significance requirements
path1_all = [WPM.pathway(BPM.path1idx)' WPM.pathway(BPM.path1idx)']';
path2_all = [WPM.pathway(BPM.path2idx)' WPM.pathway(BPM.path2idx)']';
ind_bpm = find(ismember(path1_all,pbody_pathways)|ismember(path2_all,pbody_pathways));
bpm_eff_all = [repmat({'Protective'},length(BPM.path1idx),1);repmat({'Risk'},length(BPM.path2idx),1)];

path_all = [WPM.pathway' WPM.pathway']';
ind_wpm = find(ismember(path_all,pbody_pathways));
wpm_eff_all = [repmat({'Protective'},length(WPM.ind),1);repmat({'Risk'},length(WPM.ind),1)];

ind_path = ind_wpm;
path_eff_all = wpm_eff_all;

output_bpm = cell2table(cell(0,7),'VariableNames',{'data','disease_model','path1','path2','eff_bpm','bpm_pv_discovery','bpm_ranksum_discovery'});
output_wpm = cell2table(cell(0,6),'VariableNames',{'data','disease_model','path','eff_wpm','wpm_pv_discovery','wpm_ranksum_discovery'});
output_path = cell2table(cell(0,6),'VariableNames',{'data','disease_model','path','eff_path','path_pv_discovery','path_ranksum_discovery'});

for j=1:length(dirs)
     cd(sprintf('/project/csbio/wwang/BridGE/project_PD_%s_mhygeSSI',dirs{j}))
     for i = 1:length(model)
          load(sprintf('results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',model{i}))
          ind = find(bpm_pv(ind_bpm)<=0.05);
     
          if isempty(ind)==0
               disease_model = repmat(model(i),length(ind),1);
               data = repmat(dirs(j),length(ind),1);
               path1 = path1_all(ind_bpm(ind));
               path2 = path2_all(ind_bpm(ind));
               eff_bpm = bpm_eff_all(ind_bpm(ind));
               bpm_pv_discovery = bpm_pv(ind_bpm(ind))';
               bpm_ranksum_discovery = bpm_ranksum(ind_bpm(ind))';
               output_tmp = table(data,disease_model,path1,path2,eff_bpm,bpm_pv_discovery,bpm_ranksum_discovery);
               output_bpm = [output_bpm;output_tmp];
               clear path1 path2 eff_bpm bpm_pv_discovery bpm_ranksum_discovery output_tmp
          end

          ind = find(wpm_pv(ind_wpm)<=0.05);
          if isempty(ind)==0
               disease_model = repmat(model(i),length(ind),1);
               data = repmat(dirs(j),length(ind),1);
               path = path_all(ind_wpm(ind));
               eff_wpm = wpm_eff_all(ind_wpm(ind));
               wpm_pv_discovery  = wpm_pv(ind_wpm(ind))';
               wpm_ranksum_discovery = wpm_ranksum(ind_wpm(ind))';
               output_tmp = table(data,disease_model,path,eff_wpm,wpm_pv_discovery,wpm_ranksum_discovery);
               output_wpm = [output_wpm;output_tmp];
               clear path eff_wpm wpm_pv_discovery wpm_ranksum_discovery output_tmp
          end

          ind = find(path_pv(ind_path)<=0.05);
          if isempty(ind)==0
               disease_model = repmat(model(i),length(ind),1);
               data = repmat(dirs(j),length(ind),1);
               path = path_all(ind_path(ind));
               eff_path = path_eff_all(ind_path(ind));
               path_pv_discovery = path_pv(ind_path(ind))';
               path_ranksum_discovery = path_ranksum(ind_path(ind))';
               output_tmp = table(data,disease_model,path,eff_path,path_pv_discovery,path_ranksum_discovery);
               output_path = [output_path;output_tmp];
               clear path eff_path path_pv_discovery path_ranksum_discovery output_tmp
          end
     end
end
cd('/project/csbio/wwang/BridGE/results_collection/PD')
save quickcheck_pbody_PD output_bpm output_wpm output_path

if isempty(output_bpm)==0
     writetable(output_bpm,'quickcheck_pbody_PD_bpm.xls')
end

if isempty(output_wpm)==0
     writetable(output_wpm,'quickcheck_pbody_PD_wpm.xls')
end

if isempty(output_path)==0
     writetable(output_path,'quickcheck_pbody_PD_path.xls')
end



