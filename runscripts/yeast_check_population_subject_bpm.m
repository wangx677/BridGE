% this script is used to visualize population distribution in the discovered BPMs
diseasemodel = {'RR','DD','RD','combined'};

fdrcutoff = 0.25;
pcutoff = 0.001;
BPMindfile = 'BPMind.mat';
snppathwayfile = 'snp_pathway_min5_max300.mat';

groups = readtable('/project/csbio/wwang/yeast_1011strains/maf005_groups.txt','FileType','text','ReadVariableNames',false);
groups.Var2(find(ismember(groups.Var2,''))) = {'other'};

load('SNPdataAR.mat');

for i=1:length(SNPdata.fid)
     yeast{i} = groups.Var2{find(ismember(groups.Var1,SNPdata.fid{i}))};
end
yeast = yeast';
yeast = table(yeast);
netcut = 1; % only use hygeSSI>1

for i = 1:length(diseasemodel)
     ssmFile = sprintf('ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',diseasemodel{i});
     delete(sprintf('subject_BPM_%s.xls',ssmFile))
     resultFile = sprintf('results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',diseasemodel{i});
     [subject_BPM subject_BPM_protective subject_BPM_risk summary] = summarize_bpm(ssmFile,BPMindfile,resultFile,diseasemodel{i},snppathwayfile,fdrcutoff,pcutoff,netcut);

     if isempty(subject_BPM_protective)==0
          bpm_number = sum(table2array(subject_BPM_protective(:,2:end))>0,2);
          bpm_number = table(bpm_number);
          output = [yeast subject_BPM_protective bpm_number];
          writetable(output,sprintf('subject_BPM_%s.xls',ssmFile),'Sheet',1)
     end

     if isempty(subject_BPM_risk)==0
          bpm_number = sum(table2array(subject_BPM_risk(:,2:end))>0,2);
          bpm_number = table(bpm_number);
          output = [yeast subject_BPM_risk bpm_number];
          writetable(output,sprintf('subject_BPM_%s.xls',ssmFile),'Sheet',2)
     end
end
