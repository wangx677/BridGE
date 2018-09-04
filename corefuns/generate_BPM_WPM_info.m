function generate_BPM_WPM_info(projectdir,snppathwayfile,snpgenemappingfile)

% this function is used to generate information for discovered BPMs and WPMs.
% it will one file for each BPM or WPM.

% setup directory for output files
projectdir=sprintf('%s/%s',getenv('BRIDGEPATH'),projectdir);
outputdir=sprintf('%s/BPM_WPM_info',projectdir);
if exist(outputdir,'dir')~=7
     mkdir(outputdir)
end

cd(projectdir)

models = {'RR','DD','RD','combined'};
for i=1:length(models)
     ssmfile = sprintf('ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',models{i});
     bpmindfile = 'BPMind.mat';

     try
          data = readtable(sprintf('output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',models{i}),'Sheet',3);
     catch
          data = [];
     end

     if isempty(data)~=1
          if nnz(data.fdrBPM<=0.25)>0
               data = data(find(data.fdrBPM<=0.25),:);
               for j=1:size(data,1)
                    outputfile=sprintf('%s/BPM%s_%s_%s_fdr%.2f_pv%.2e_rs%.2f',outputdir,num2str(j),models{i},data.eff_bpm{j},data.fdrBPM(j),data.bpm_pv_discovery(j),data.bpm_ranksum_discovery(j));
                    get_interaction_pair(data.path1{j},data.path2{j},data.eff_bpm{j},ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,outputfile)
               end
          end
     end

     try          
          data = readtable(sprintf('output_results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat.xls',models{i}),'Sheet',4);
     catch
          data = [];
     end

     if isempty(data)~=1
          if nnz(data.fdrWPM<=0.25)>0
               data = data(find(data.fdrWPM<=0.25),:);
               for j=1:size(data,1)
                    outputfile=sprintf('%s/WPM%s_%s_%s_fdr%.2f_pv%.2e_rs%.2f',outputdir,num2str(j),models{i},data.eff_wpm{j},data.fdrWPM(j),data.wpm_pv_discovery(j),data.wpm_ranksum_discovery(j));
                    get_interaction_pair(data.path{j},data.path{j},data.eff_wpm{j},ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,outputfile)
               end
          end
     end
end
