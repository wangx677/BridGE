function density2pvalue(ssmFile,samplePerms)

% Inputs
% ssmFile: SNP-SNP interaction file
% samplePerms: number of sample permutation
%
% Outputs:
% bpm_pv, wpm_pv: permutation p-values

for i=0:samplePerms
     load(sprintf('genstats_%s_R%s.mat',ssmFile,num2str(i)))
     for tt=1:2
          bpm_density_all{tt}(i+1,:) = BPM_density{tt};
          wpm_density_all{tt}(i+1,:) = WPM_density{tt};
          
          bpm_local_all{tt}(i+1,:) = bpm_local{tt};
          wpm_local_all{tt}(i+1,:) = wpm_local{tt};
     end
end

load BPMind.mat
     
for i=1:samplePerms+1
     for tt=1:2
          ind = setdiff(1:samplePerms+1,i);
          bpm_pv_density{tt} = (sum(bpm_density_all{tt}(ind,:)>=bpm_density_all{tt}(i,:))+1)/(samplePerms+1);
          wpm_pv_density{tt} = (sum(wpm_density_all{tt}(ind,:)>=wpm_density_all{tt}(i,:))+1)/(samplePerms+1);

          bpm_pv_density_and_ranksum{tt} = (sum(bpm_density_all{tt}(ind,:)>=bpm_density_all{tt}(i,:) & bpm_local_all{tt}(ind,:)>=bpm_local_all{tt}(i,:))+1)/(samplePerms+1);
          wpm_pv_density_and_ranksum{tt} = (sum(wpm_density_all{tt}(ind,:)>=wpm_density_all{tt}(i,:) & wpm_local_all{tt}(ind,:)>=wpm_local_all{tt}(i,:))+1)/(samplePerms+1);


          % assign 1 to pv when the density is nan
          bpm_pv_density{tt}(isnan(bpm_density_all{tt}(i,:))) = 1;
          wpm_pv_density{tt}(isnan(wpm_density_all{tt}(i,:))) = 1;

          bpm_pv_density_and_ranksum{tt}(isnan(bpm_density_all{tt}(i,:))) = 1;
          wpm_pv_density_and_ranksum{tt}(isnan(wpm_density_all{tt}(i,:))) = 1;

          % if BPM size is less than 5x5, set pv to be 1
          bpm_pv_density{tt}(find(BPM.ind1size<5 | BPM.ind2size<5)) = 1;
          bpm_pv_density_and_ranksum{tt}(find(BPM.ind1size<5 | BPM.ind2size<5)) = 1;

          % if ranksum test p-value is >0.05, set bpm/wpm_pv_density to be 1
          bpm_pv_density{tt}(find(bpm_local_all{tt}(i,:)<-log10(0.05))) = 1;
          wpm_pv_density{tt}(find(wpm_local_all{tt}(i,:)<-log10(0.05))) = 1;
     end

     save(sprintf('genstats_%s_R%s.mat',ssmFile,num2str(i-1)),'bpm_pv*','wpm_pv*','-append')
end
