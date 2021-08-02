function [number_GI_risk, number_GI_validated_risk, number_GI_protective, number_GI_validated_protective] = summarize_result_gene(projectdirs,ssmFile,R,fdrcutoff,pvcutoff,outputdir,validcutoff,tail,goi)

% projectdirs: bridge directories withou summary results files
% GI: hygeSSI or mhygeSSI
% fdrcutoff: fdr cutoff value, fdrcutoff=1 means using pvcutoff
% pvcutoff: p-value cutoff value, pvcutoff=1 means using fdrcutoff
% outputdir: output files directory
% validcutoff: number of datasets that can validate the discovery
% tail: one-tail or two-tail test
% goi: gene of interest. If 'goi' exist, only count goi genes assoicated interactions

union_risk_pair = array2table(zeros(0,2));
union_risk_pair.Properties.VariableNames = {'gene1','gene2'};

union_protective_pair = array2table(zeros(0,2));
union_protective_pair.Properties.VariableNames = {'gene1','gene2'};

for i=1:length(projectdirs)
     load(sprintf('/project/csbio/wwang/BridGE/%s/summary_results_%s_R%s.mat',projectdirs{i},ssmFile,num2str(R)))
     result_protective_full{i} = result_protective_all;
     result_risk_full{i} = result_risk_all;

     union_protective_pair = [union_protective_pair;result_protective_all(:,1:2)];
     union_risk_pair = [union_risk_pair;result_risk_all(:,1:2)];
end

union_protective_pair = unique(union_protective_pair);
union_risk_pair = unique(union_risk_pair);

union_result_protective = array2table(nan(size(union_protective_pair,1),length(projectdirs)*2));
union_result_risk = array2table(nan(size(union_risk_pair,1),length(projectdirs)*2));

for i=1:length(projectdirs)
     [tmp, ind] = ismember(result_protective_full{i}(:,1:2),union_protective_pair(:,1:2),'rows');
     if strcmp(tail,'one')
          union_result_protective(ind,2*i-1:2*i) = result_protective_full{i}(:,[6 7]);
     elseif strcmp(tail,'two')
          union_result_protective(ind,2*i-1:2*i) = result_protective_full{i}(:,[6 8]);
     end

     [tmp, ind] = ismember(result_risk_full{i}(:,1:2),union_risk_pair(:,1:2),'rows');
     if strcmp(tail,'one')
          union_result_risk(ind,2*i-1:2*i) = result_risk_full{i}(:,[6 7]); 
     elseif strcmp(tail,'two')
          union_result_risk(ind,2*i-1:2*i) = result_risk_full{i}(:,[6 8]);
     end 
end

for i=1:(size(union_result_protective,2)/2)
     varnames{2*i-1} = sprintf('pv_data%s',num2str(i));
     varnames{2*i} = sprintf('fdr_data%s',num2str(i));
end

union_result_protective.Properties.VariableNames = varnames;
union_result_risk.Properties.VariableNames = varnames;

effect = repmat({'protective'},size(union_result_protective,1),1);
effect = table(effect);
union_result_protective = [union_protective_pair effect union_result_protective];

effect = repmat({'risk'},size(union_result_risk,1),1);
effect = table(effect);
union_result_risk = [union_risk_pair effect union_result_risk];

% if "goi" exist
if exist('goi')
     ind = find(ismember(union_result_protective.gene1,goi)|ismember(union_result_protective.gene2,goi));
     union_result_protective = union_result_protective(ind,:);

     ind = find(ismember(union_result_risk.gene1,goi)|ismember(union_result_risk.gene2,goi));
     union_result_risk = union_result_risk(ind,:);
end

% find number of validated interactions (validated in 1, 2,..,n-1 dataset)
clear ind

for i=1:length(projectdirs);
     if pvcutoff==1
          ind1_tmp = find(table2array(union_result_protective(:,3+2*i))<=fdrcutoff);
     elseif fdrcutoff==1
          ind1_tmp = find(table2array(union_result_protective(:,3+2*i-1))<=pvcutoff);
     end
     I = setdiff(1:length(projectdirs),i);
     ind2_tmp = find(sum(table2array(union_result_protective(:,3+I*2-1))<=0.05,2)>=validcutoff);
     ind_protective{i} = intersect(ind1_tmp,ind2_tmp);

     if pvcutoff==1
          ind1_tmp = find(table2array(union_result_risk(:,3+2*i))<=fdrcutoff);
     elseif fdrcutoff==1
          ind1_tmp = find(table2array(union_result_risk(:,3+2*i-1))<=pvcutoff);
     end
     I = setdiff(1:length(projectdirs),i);
     ind2_tmp = find(sum(table2array(union_result_risk(:,3+I*2-1))<=0.05,2)>=validcutoff);
     ind_risk{i} = intersect(ind1_tmp,ind2_tmp);

end

if sum(cellfun(@(x)isempty(x),ind_protective))<=0
     number_GI_validated_protective = length(unique(cell2mat(ind_protective')));
else
     ind = find(cellfun(@(x)isempty(x),ind_protective)~=1);
     if length(ind)>0
           number_GI_validated_protective = length(unique(cell2mat(ind_protective(ind)')));
     else
          number_GI_validated_protective = 0;
     end
end

if sum(cellfun(@(x)isempty(x),ind_risk))<=0
     number_GI_validated_risk = length(unique(cell2mat(ind_risk')));
else
     ind = find(cellfun(@(x)isempty(x),ind_risk)~=1);
     if length(ind)>0
           number_GI_validated_risk = length(unique(cell2mat(ind_risk(ind)')));
     else
          number_GI_validated_risk = 0;
     end
end

% only keep gene pairs that have p<=0.05 in at least one cohort
I = 1:length(projectdirs);
ind = find(min(table2array(union_result_protective(:,3+I*2-1)),[],2)<=0.05);
union_result_protective = union_result_protective(ind,:);

% get number of GI's pass FDR cutoff or pv cutoff
if pvcutoff==1
     number_GI_protective = length(find(min(table2array(union_result_protective(:,3+I*2)),[],2)<=fdrcutoff));
elseif fdrcutoff==1
      number_GI_protective = length(find(min(table2array(union_result_protective(:,3+I*2-1)),[],2)<=pvcutoff));
end

% summaryize  minimum FDR, corresponding p-value, and numbe of validation for each discovery
if size(union_result_protective,1)>0
     [minFDR  idx]= min(table2array(union_result_protective(:,3+I*2)),[],2);
     pv = table2array(union_result_protective(:,3+I*2-1));
     pv1 = arrayfun(@(x)pv(x,idx(x)),1:size(pv,1));
     pv = reshape(pv1,length(pv1),1);
     number_validation = sum(table2array(union_result_protective(:,3+I*2-1))<=0.05,2)-1;
     number_validation(find(minFDR>fdrcutoff)) = 0;
     minFDR = table(minFDR);
     pv = table(pv);
     number_validation = table(number_validation);
     union_result_protective = [union_result_protective minFDR pv number_validation];
     [tmp idx] = sortrows([union_result_protective.number_validation union_result_protective.minFDR],[-1 2]);
     union_result_protective = union_result_protective(idx,:);
end

% repeat above for risk interactions
ind = find(min(table2array(union_result_risk(:,3+I*2-1)),[],2)<=0.05);
union_result_risk = union_result_risk(ind,:);

if pvcutoff==1
     number_GI_risk = length(find(min(table2array(union_result_risk(:,3+I*2)),[],2)<=fdrcutoff));
elseif fdrcutoff==1
     number_GI_risk = length(find(min(table2array(union_result_risk(:,3+I*2-1)),[],2)<=pvcutoff));
end

if size(union_result_risk,1) > 0
     [minFDR  idx]= min(table2array(union_result_risk(:,3+I*2)),[],2);
     pv = table2array(union_result_risk(:,3+I*2-1));
     pv1 = arrayfun(@(x)pv(x,idx(x)),1:size(pv,1));
     pv = reshape(pv1,length(pv1),1);
     number_validation = sum(table2array(union_result_risk(:,3+I*2-1))<=0.05,2)-1;
     number_validation(find(minFDR>fdrcutoff)) = 0;
     minFDR = table(minFDR);
     pv = table(pv);
     number_validation = table(number_validation);
     union_result_risk = [union_result_risk minFDR pv number_validation];
     [tmp idx] = sortrows([union_result_risk.number_validation union_result_risk.minFDR],[-1 2]);
     union_result_risk = union_result_risk(idx,:);
end

if R==0 
     cd(outputdir)
     
     if exist('goi')
          if pvcutoff==1 & fdrcutoff==0.25
               writetable(union_result_protective,sprintf('results_GOI_%s_%s_FDR%s_validatedin%s_%stailed.xls',goi,ssmFile,num2str(fdrcutoff),num2str(validcutoff),tail),'Sheet',1)
               writetable(union_result_risk,sprintf('results_GOI_%s_%s_FDR%s_validatedin%s_%stailed.xls',goi,ssmFile,num2str(fdrcutoff),num2str(validcutoff),tail),'Sheet',2)
          elseif fdrcutoff==1 & pvcutoff==0.05
               writetable(union_result_protective,sprintf('results_GOI_%s_%s_pv%s_validatedin%s_%stailed.xls',goi,ssmFile,num2str(pvcutoff),num2str(validcutoff),tail),'Sheet',1)
               writetable(union_result_risk,sprintf('results_GOI_%s_%s_pv%s_validatedin%s_%stailed.xls',goi,ssmFile,num2str(pvcutoff),num2str(validcutoff),tail),'Sheet',2)
          end
     else
          if pvcutoff==1 & fdrcutoff==0.25
               writetable(union_result_protective,sprintf('results_GOI_%s_FDR%s_validatedin%s_%stailed.xls',ssmFile,num2str(fdrcutoff),num2str(validcutoff),tail),'Sheet',1)
               writetable(union_result_risk,sprintf('results_GOI_%s_FDR%s_validatedin%s_%stailed.xls',ssmFile,num2str(fdrcutoff),num2str(validcutoff),tail),'Sheet',2)
          elseif fdrcutoff==1 & pvcutoff==0.05
                writetable(union_result_protective,sprintf('results_GOI_%s_pv%s_validatedin%s_%stailed.xls',ssmFile,num2str(pvcutoff),num2str(validcutoff),tail),'Sheet',1)
               writetable(union_result_risk,sprintf('results_GOI_%s_pv%s_validatedin%s_%stailed.xls',ssmFile,num2str(pvcutoff),num2str(validcutoff),tail),'Sheet',2)
          end
     end
end
 
