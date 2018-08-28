function update_wpmdensity(model,N,snpPerms)

% add wpm_density to result files because there was a typo in genstats_zscore.m 
load BPMind.mat
BPMsize =[BPM.size BPM.size];
WPMsize = [WPM.size WPM.size];
WPMind = WPM.ind;
BPMpath1size = [BPM.ind1size BPM.ind1size];
BPMpath2size = [BPM.ind2size BPM.ind2size];

clear BPM WPM
for i=0:N
     file=sprintf('genstats_zscore_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R%s_snpPerm%s.mat',model,num2str(i),num2str(snpPerms));
     load(file)

     if exist('wpm_density','var')==0
          load(sprintf('ssM_hygeSSI_alpha10.05_alpha20.05_%s_R%s.mat',model,num2str(i)))
          if min(size(ssM{1},1),size(ssM{1},2))==1
              for tt=1:2
                    ssM{tt} = squareform(ssM{tt});
               end
          end
          for tt=1:2
               wpm_density{tt} = cellfun(@(x)sum2(ssM{tt}(x,x)),WPMind);
               wpm_density{tt} = wpm_density{tt}./WPMsize(1:length(WPMsize)/2);
          end
          save(file,'wpm_density','-append')
          clear ssM
     end
     clear bpm_pv_snp wpm_pv_snp path_pv_snp bpm_density path_density wpm_density ind2keep_bpm maxidx minidx bpm_zscore_snp wpm_zscore_snp path_zscore_snp
end

