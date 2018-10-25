function summary(bridgePATH,dirlist,outputfile)

N = length(dirlist);
output = [];

for i=1:N

     pattern = fullfile(sprintf('%s/%s',bridgePATH,dirlist{i}), 'result*.mat');
     files = dir(pattern);

     for k = 1:length(files)
          fileName = fullfile(sprintf('%s/%s',bridgePATH,dirlist{i}), sprintf('output_%s',files(k).name));
          phenotype = strsplit(dirlist{i},'_');
          phenotype = strjoin(phenotype(2:end),'_');
          phenotype = repmat({phenotype},3,1);
          diseaseModel = strsplit(files(k).name,'_');
          diseaseModel = diseaseModel(end-1);
          diseaseModel = repmat(diseaseModel,3,1);
          GItype = {'BPM';'WPM';'PATH'};

          if exist(fileName,'file')==2
               data = load(fileName);
               outputtmp = [table(phenotype,diseaseModel,GItype) data.output_discovery_summary];
               outputtmp.Properties.RowNames={};
          else
               outputtmp = array2table(zeros(3,size(output,2)-3),'VariableNames',{'fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40'}); 
               outputtmp = [table(phenotype,diseaseModel,GItype) outputtmp];
          end
          output = [output;outputtmp];
     end
end
writetable(output,outputfile)     
