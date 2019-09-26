function collectresults(resultfile,fdrcut,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,validationfile)
% resultfile='/project/csbio/wwang/BridGE_nobin/project_PD_Simon/results_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat';
% validationfile='/project/csbio/wwang/BridGE_nobin/project_PD_NGRC/results_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat';

diseasemodel = strsplit(resultfile,'_');
diseasemodel = diseasemodel{end-1};

% load resultfile and set FDR to be 1 if no FDR exist
load(resultfile)

if exist('fdrBPM2','var')
     fdrBPM = fdrBPM2;
else
     fdrBPM = ones(1,length(bpm_pv))';
end

if exist('fdrWPM2','var') 
     fdrWPM = fdrWPM2;
else
     if exist('wpm_pv','var')
          fdrWPM = ones(1,length(wpm_pv));
     end
end

if exist('fdrPATH2','var')
     fdrPATH = fdrPATH2;
else
     if exist('path_pv','var')
          fdrPATH = ones(1,length(path_pv));
     end
end

% find BPM/WPM/PATH that passed the FDR cutoff
if exist('fdrBPM','var')
     ind_bpm = find(fdrBPM<=fdrcut);
else
     ind_bpm = [];     
end

if exist('fdrWPM','var')     
     ind_wpm = find(fdrWPM<=fdrcut);
else
     ind_wpm = [];
end

if exist('fdrPATH','var')     
     ind_path = find(fdrPATH<=fdrcut);
else
     ind_path = [];
end

% if there exist BPM/WPM/PATH passed the cutoff
% collect corresponding stats

if length(ind_bpm>0)
     fdrBPM = fdrBPM(ind_bpm)';
     bpm_pv_discovery = bpm_pv(ind_bpm)';
     bpm_ranksum_discovery = bpm_ranksum(ind_bpm)';
     
     % retreive relevant information about effect size
     netcut = 1;
     pcutoff = 1;
     [~, ~, ~, effectsize_bpm] = summarize_bpm(ssmfile,bpmindfile,resultfile,diseasemodel,snppathwayfile,fdrcut,pcutoff,netcut);
end

if (length(ind_wpm)>0)
     fdrWPM = fdrWPM(ind_wpm)';
     wpm_pv_discovery = wpm_pv(ind_wpm)';
     wpm_ranksum_discovery = wpm_ranksum(ind_wpm)';
     
     % retreive relevant information about effect size
     netcut = 1;
     pcutoff = 1;
     [~, ~, ~, effectsize_wpm] = summarize_wpm(ssmfile,bpmindfile,resultfile,diseasemodel,snppathwayfile,fdrcut,pcutoff,netcut);
end

if (length(ind_path)>0)     
     fdrPATH = fdrPATH(ind_path)';
     path_pv_discovery = path_pv(ind_path)';
     path_ranksum_discovery = path_ranksum(ind_path)';
end

clear fdrBPM2 fdrWPM2 fdrPATH2 bpm_pv wpm_pv path_pv bpm_ranksum wpm_ranksum path_ranksum

% if validationfile exist, collect validation stats
if exist('validationfile','var')==1
     load(validationfile)
     if length(ind_bpm>0)
          bpm_ranksum_valid = bpm_ranksum(ind_bpm)';
          bpm_pv_valid = bpm_pv(ind_bpm)';
     end

     if (length(ind_wpm)>0)
          wpm_ranksum_valid = wpm_ranksum(ind_wpm)';
          wpm_pv_valid = wpm_pv(ind_wpm)';
     end

     if (length(ind_path)>0)
          path_ranksum_valid = path_ranksum(ind_path)';
          path_pv_valid = path_pv(ind_path)';
     end
end
clear fdrBPM2 fdrWPM2 fdrPATH2 bpm_pv wpm_pv path_pv bpm_ranksum wpm_ranksum path_ranksum

% retreive relevant information for driver SNPs and corresponding genes
if (length(ind_bpm)>0 | length(ind_wpm)>0 | length(ind_path)>0)
     load(bpmindfile);
     load(snppathwayfile);
     genesetfile = snpset.genesets;
     genesetfile = sprintf('../refdata/%s',genesetfile);

     if length(ind_bpm)>0
          if isfield(BPM,'path1')==0 
               path1 = snpset.pathwaynames(BPM.path1idx);
               path2 = snpset.pathwaynames(BPM.path2idx);
          else
               path1 = BPM.path1;
               path2 = BPM.path2;
          end
          path1 = repmat(reshape(path1,length(path1),1),2,1);
          path2 = repmat(reshape(path2,length(path2),1),2,1);

          path1 = path1(ind_bpm);
          path2 = path2(ind_bpm);

          bpm_size = repmat(reshape(BPM.size,length(BPM.size),1),2,1);
          bpm_size = bpm_size(ind_bpm);

          eff_bpm = cell(size(path1));
          for i=1:length(ind_bpm)
               if ind_bpm(i)<=length(BPM.size)
                    eff_bpm{i} = 'protective';
               else
                    eff_bpm{i} = 'risk';
               end
          end

          eff_bpm = reshape(eff_bpm,length(eff_bpm),1);
                    
          for i = 1:length(ind_bpm)
               [output_path1_snp, output_path2_snp] = get_interaction_pair(path1{i},path2{i},eff_bpm{i},ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile);
               idx = find(output_path1_snp.gi_fold>1);
               
               if length(idx)>0
                    for j=1:min(20,length(idx))
                         bpm_snp1_tmp{j} = strjoin({output_path1_snp.snps{j},output_path1_snp.genes{j},sprintf('%.2f',output_path1_snp.gi_fold(j))},'_');
                    end

                    bpm_path1_drivers{i} = strjoin(bpm_snp1_tmp,';');
               else
                     bpm_path1_drivers{i} = nan;
               end

               idx = find(output_path2_snp.gi_fold>1);
               if length(idx)>0
                    for j=1:min(20,length(idx))
                         bpm_snp2_tmp{j} = strjoin({output_path2_snp.snps{j},output_path2_snp.genes{j},sprintf('%.2f',output_path2_snp.gi_fold(j))},'_');
                    end
                    
                    bpm_path2_drivers{i} = strjoin(bpm_snp2_tmp,';');
               else
                    bpm_path2_drivers{i} = nan;
               end     
  
               clear output_path1_snp output_path2_snp bpm_snp1_tmp bpm_snp2_tmp
          end

          bpm_path1_drivers = reshape(bpm_path1_drivers,length(bpm_path1_drivers),1);
          bpm_path2_drivers = reshape(bpm_path2_drivers,length(bpm_path2_drivers),1);
     end
     
     if length(ind_wpm)>0     
          pathway = repmat(reshape(snpset.pathwaynames,length(snpset.pathwaynames),1),2,1);
          path_wpm = pathway(ind_wpm);
          
          wpm_size = repmat(reshape(WPM.size,length(WPM.size),1),2,1);
          wpm_size = wpm_size(ind_wpm);

          eff_wpm = cell(size(wpm_size));
          for i=1:length(ind_wpm)
               if ind_wpm(i)<=length(WPM.size)
                   eff_wpm{i} = 'protective';
               else
                    eff_wpm{i} = 'risk';      
               end
          end

          eff_wpm = reshape(eff_wpm,length(eff_wpm),1);

          for i=1:length(ind_wpm)
               [output_path_snp] = get_interaction_pair(path_wpm{i},path_wpm{i},eff_wpm{i},ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile);
               idx = find(output_path_snp.gi_fold>1);
               
               if length(idx)>0
                    for j=1:min(20,length(idx))
                         wpm_snp_tmp{j} = strjoin({output_path_snp.snps{j},output_path_snp.genes{j},sprintf('%.2f',output_path_snp.gi_fold(j))},'_');
                    end
               
                    wpm_path_drivers{i} = strjoin(wpm_snp_tmp,';');
               else
                    wpm_path_drivers{i} = nan;
               end

               clear output_path_snp wpm_snp_tmp
          end
          
          wpm_path_drivers = reshape(wpm_path_drivers,length(wpm_path_drivers),1);
     end

     if length(ind_path)>0
          pathway = repmat(reshape(snpset.pathwaynames,length(snpset.pathwaynames),1),2,1);
          path_path = pathway(ind_path);

          path_size = repmat(reshape(WPM.indsize,length(WPM.indsize),1),2,1);
          path_size = path_size(ind_path);

          eff_path = cell(size(path_size));
          for i=1:length(ind_path)
               if ind_path(i)<=length(WPM.size)
                    eff_path{i} = 'protective';
               else
                    eff_path{i} = 'risk';
               end
          end

          eff_path = reshape(eff_path,length(eff_path),1);
     end
end

% prepare output table
if exist('validationfile','var')==1
     if length(ind_bpm)>0
          fdrBPM = round(fdrBPM*100)/100;
          bpm_ranksum_discovery = round(bpm_ranksum_discovery*100)/100;
          bpm_ranksum_valid = round(bpm_ranksum_valid*100)/100;
          bpm_eff_case = round(effectsize_bpm.cases*100)/100;
          bpm_eff_control = round(effectsize_bpm.controls*100)/100;
          bpm_eff_size_OR = round(effectsize_bpm.oddsratio*100)/100; 
          output_bpm_table = table(path1,path2,fdrBPM,eff_bpm,bpm_size,bpm_pv_discovery,bpm_ranksum_discovery,bpm_eff_case,bpm_eff_control,bpm_eff_size_OR,bpm_pv_valid,bpm_ranksum_valid,bpm_path1_drivers,bpm_path2_drivers);
     end

     if length(ind_wpm)>0
          path = path_wpm;
          fdrWPM = round(fdrWPM*100)/100;
          wpm_ranksum_discovery = round(wpm_ranksum_discovery*100)/100;
          wpm_ranksum_valid = round(wpm_ranksum_valid*100)/100;
          wpm_eff_case = round(effectsize_wpm.cases*100)/100;
          wpm_eff_control = round(effectsize_wpm.controls*100)/100;
          wpm_eff_size_OR = round(effectsize_wpm.oddsratio*100)/100;
          output_wpm_table = table(path,fdrWPM,eff_wpm,wpm_size,wpm_pv_discovery,wpm_ranksum_discovery,wpm_eff_case,wpm_eff_control,wpm_eff_size_OR,wpm_pv_valid,wpm_ranksum_valid,wpm_path_drivers);
     end

     if length(ind_path)>0
          path = path_path;
          fdrPATH = round(fdrPATH*100)/100;
          path_ranksum_discovery = round(path_ranksum_discovery*100)/100;
          path_ranksum_valid = round(path_ranksum_valid*100)/100;
          output_path_table = table(path,fdrPATH,eff_path,path_size,path_pv_discovery,path_ranksum_discovery,path_pv_valid,path_ranksum_valid);
     end

else
     if length(ind_bpm)>0
          fdrBPM = round(fdrBPM*100)/100;
          bpm_ranksum_discovery = round(bpm_ranksum_discovery*100)/100;
          bpm_eff_case = round(effectsize_bpm.cases*100)/100;
          bpm_eff_control = round(effectsize_bpm.controls*100)/100;
          bpm_eff_size_OR = round(effectsize_bpm.oddsratio*100)/100;
          output_bpm_table = table(path1,path2,fdrBPM,eff_bpm,bpm_size,bpm_pv_discovery,bpm_ranksum_discovery,bpm_eff_case,bpm_eff_control,bpm_eff_size_OR,bpm_path1_drivers,bpm_path2_drivers);
     end

     if length(ind_wpm)>0
          path = path_wpm;
          fdrWPM = round(fdrWPM*100)/100;
          wpm_ranksum_discovery = round(wpm_ranksum_discovery*100)/100;
          wpm_eff_case = round(effectsize_wpm.cases*100)/100;
          wpm_eff_control = round(effectsize_wpm.controls*100)/100;
          wpm_eff_size_OR = round(effectsize_wpm.oddsratio*100)/100;
          output_wpm_table = table(path,fdrWPM,eff_wpm,wpm_size,wpm_pv_discovery,wpm_ranksum_discovery,wpm_eff_case,wpm_eff_control,wpm_eff_size_OR,wpm_path_drivers);
     end

     if length(ind_path)>0
          path = path_path;
          fdrPATH = round(fdrPATH*100)/100;
          path_ranksum_discovery = round(path_ranksum_discovery*100)/100;
          output_path_table = table(path,fdrPATH,eff_path,path_size,path_pv_discovery,path_ranksum_discovery);
     end
end

% discovery summary: count number of signficance discoveries at given cutoffs
output_discovery_summary(1,:) = [min(fdrBPM) nnz(fdrBPM<=0.05) nnz(fdrBPM<=0.1) nnz(fdrBPM<=0.15) nnz(fdrBPM<=0.2) nnz(fdrBPM<=0.25) nnz(fdrBPM<=0.3) nnz(fdrBPM<=0.35) nnz(fdrBPM<=0.4)];    
if exist('fdrWPM','var')
     output_discovery_summary(2,:) = [min(fdrWPM) nnz(fdrWPM<=0.05) nnz(fdrWPM<=0.1) nnz(fdrWPM<=0.15) nnz(fdrWPM<=0.2) nnz(fdrWPM<=0.25) nnz(fdrWPM<=0.3) nnz(fdrWPM<=0.35) nnz(fdrWPM<=0.4)];
     output_discovery_summary(3,:) = [min(fdrPATH) nnz(fdrPATH<=0.05) nnz(fdrPATH<=0.1) nnz(fdrPATH<=0.15) nnz(fdrPATH<=0.2) nnz(fdrPATH<=0.25) nnz(fdrPATH<=0.3) nnz(fdrPATH<=0.35) nnz(fdrPATH<=0.4)];
end

% non-redundant discovery summary: count number of non-redundant signficance discoveries at given cutoffs
[BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group_tmp,WPM_group_tmp,PATH_group_tmp] = check_BPM_WPM_redundancy(fdrBPM,fdrWPM,fdrPATH,bpmindfile,0.4);
output_noRD_discovery_summary(1,:) = [min(fdrBPM) BPM_nosig_noRD];

if exist('fdrWPM','var')
     output_noRD_discovery_summary(2,:) = [min(fdrWPM) WPM_nosig_noRD];
     output_noRD_discovery_summary(3,:) = [min(fdrPATH) PATH_nosig_noRD];
end

% 
if exist('validationfile','var')==1
     if exist('bpm_ranksum_valid','var')==1
          output_validation_summary(1,:) = [nnz(fdrBPM<=0.05 & bpm_ranksum_valid>=-log10(0.05)) nnz(fdrBPM<=0.1 & bpm_ranksum_valid>=-log10(0.05)) ...
                         nnz(fdrBPM<=0.15 & bpm_ranksum_valid>=-log10(0.05)) nnz(fdrBPM<=0.2 & bpm_ranksum_valid>=-log10(0.05)) nnz(fdrBPM<=0.25 & bpm_ranksum_valid>=-log10(0.05)) ...
                          nnz(fdrBPM<=0.3 & bpm_ranksum_valid>=-log10(0.05)) nnz(fdrBPM<=0.35 & bpm_ranksum_valid>=-log10(0.05)) nnz(fdrBPM<=0.4 & bpm_ranksum_valid>=-log10(0.05))];
     else
          output_validation_summary(1,:) = [0 0 0 0 0 0 0 0];
     end

     if exist('wpm_ranksum_valid','var')==1
          output_validation_summary(2,:) = [nnz(fdrWPM<=0.05 & wpm_ranksum_valid>=-log10(0.05)) nnz(fdrWPM<=0.1 & wpm_ranksum_valid>=-log10(0.05)) ...
                         nnz(fdrWPM<=0.15 & wpm_ranksum_valid>=-log10(0.05)) nnz(fdrWPM<=0.2 & wpm_ranksum_valid>=-log10(0.05)) nnz(fdrWPM<=0.25 & wpm_ranksum_valid>=-log10(0.05)) ...
                         nnz(fdrWPM<=0.3 & wpm_ranksum_valid>=-log10(0.05)) nnz(fdrWPM<=0.35 & wpm_ranksum_valid>=-log10(0.05)) nnz(fdrWPM<=0.4 & wpm_ranksum_valid>=-log10(0.05))];
     else
          output_validation_summary(2,:) = [0 0 0 0 0 0 0 0];
     end

     if exist('path_ranksum_valid','var')==1
          output_validation_summary(3,:) = [nnz(fdrPATH<=0.05 & path_ranksum_valid>=-log10(0.05)) nnz(fdrPATH<=0.1 & path_ranksum_valid>=-log10(0.05)) ...
                         nnz(fdrPATH<=0.15 & path_ranksum_valid>=-log10(0.05)) nnz(fdrPATH<=0.2 & path_ranksum_valid>=-log10(0.05)) nnz(fdrPATH<=0.25 & path_ranksum_valid>=-log10(0.05)) ...
                         nnz(fdrPATH<=0.3 & path_ranksum_valid>=-log10(0.05)) nnz(fdrPATH<=0.35 & path_ranksum_valid>=-log10(0.05)) nnz(fdrPATH<=0.4 & path_ranksum_valid>=-log10(0.05))];
     else
          output_validation_summary(3,:) = [ 0 0 0 0 0 0 0 0];
     end

     output_validation_summary = round((output_validation_summary./output_discovery_summary)*100)/100;

     output_validation_summary = array2table(output_validation_summary,'VariableNames',{'fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40'},'RowNames',{'BPM','WPM','PATH'});
end

if exist('fdrWPM','var')
     output_discovery_summary = array2table(output_discovery_summary,'VariableNames',{'minfdr','fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40'},'RowNames',{'BPM','WPM','PATH'});

     output_noRD_discovery_summary = array2table(output_noRD_discovery_summary,'VariableNames',{'minfdr','fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40'},'RowNames',{'BPM','WPM','PATH'});
else
     output_discovery_summary = array2table(output_discovery_summary,'VariableNames',{'minfdr', 'fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40'},'RowNames',{'BPM'});
     output_noRD_discovery_summary = array2table(output_noRD_discovery_summary,'VariableNames',{'minfdr','fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40'},'RowNames',{'BPM'});
end

filename = strsplit(resultfile,'/');
filename = filename{end};


if exist('validationfile','var')==1
     save(sprintf('output_%s',filename),'output*','resultfile','validationfile');
else
     save(sprintf('output_%s',filename),'output*','resultfile');
end

filename = sprintf('output_%s.xls',filename);

if exist(filename,'file')==2
     delete(filename)
end

% sort results
if exist('bpm_pv_discovery','var')
     [tmp ind_bpm] = sortrows([fdrBPM bpm_pv_discovery bpm_ranksum_discovery],[1 2 -3]);
     output_bpm_table = output_bpm_table(ind_bpm,:);
     index = [1:size(output_bpm_table,1)]';
     index = table(index);
     output_bpm_table = [index output_bpm_table];
end

if exist('wpm_pv_discovery','var')
     [tmp ind_wpm] = sortrows([fdrWPM wpm_pv_discovery wpm_ranksum_discovery],[1 2 -3]);
     output_wpm_table = output_wpm_table(ind_wpm,:);
     index = [1:size(output_wpm_table,1)]';
     index = table(index);
     output_wpm_table = [index output_wpm_table];
end

if exist('path_pv_discovery','var')
     [tmp ind_path] = sortrows([fdrPATH path_pv_discovery path_ranksum_discovery],[1 2 -3]);
     output_path_table = output_path_table(ind_path,:);
     index = [1:size(output_path_table,1)]';
     index = table(index);
     output_path_table = [index output_path_table];
end

writetable(output_discovery_summary,filename,'Sheet',1,'WriteRowNames',true)

if exist('validationfile','var')==1
     writetable(output_validation_summary,filename,'Sheet',2,'WriteRowNames',true)
end

if length(ind_bpm)>0
     writetable(output_bpm_table,filename,'Sheet',3)
end

if length(ind_wpm)>0
     writetable(output_wpm_table,filename,'Sheet',4)
end

if length(ind_path)>0
     writetable(output_path_table,filename,'Sheet',5)
end

