function run_analysis_randomize_SNP(model,R)

load BPMind.mat
     
load(sprintf('ssM_mhygeSSI_alpha10.05_alpha20.05_%s_R0.mat',model));
for tt=1:2
     ssM{tt} = squareform(ssM{tt});
end

rng(R)
ind = randperm(size(ssM{tt},1));
for tt=1:2
     ssM{tt} = ssM{tt}(ind,ind);
end

for tt=1:2
     BPM_density{tt}(1,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM{tt}(x,y),1,[])),BPM.ind1,BPM.ind2,'uniform',0));
     WPM_density{tt}(1,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM{tt}(x,y),1,[])),WPM.ind,WPM.ind,'uniform',0));
end

for i=1:1000
     load(sprintf('ssM_mhygeSSI_alpha10.05_alpha20.05_%s_R%s.mat',model,num2str(i)));
     for tt=1:2
          ssM{tt} = squareform(ssM{tt});
          ssM{tt} = ssM{tt}(ind,ind);
     end

     for tt=1:2
          tmp_BPM{tt}(i,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM{tt}(x,y),1,[])),BPM.ind1,BPM.ind2,'uniform',0));
          tmp_WPM{tt}(i,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM{tt}(x,y),1,[])),WPM.ind,WPM.ind,'uniform',0));
     end
end


for tt=1:2
     BPM_density_pv{tt} = (arrayfun(@(x)nnz(tmp_BPM{tt}(:,x)>=BPM_density{tt}(x)),1:length(BPM_density{tt}))+1)/1000;
     WPM_density_pv{tt} = (arrayfun(@(x)nnz(tmp_WPM{tt}(:,x)>=WPM_density{tt}(x)),1:length(WPM_density{tt}))+1)/1000;
end

for tt=1:2
     BPM_density_pv{tt}(isnan(BPM_density{tt})) = 1;
     WPM_density_pv{tt}(isnan(WPM_density{tt})) = 1;
end

% compute fdr
for tt=1:2
     fdr_BPM{tt} = pvalue_rank_fdr(BPM_density_pv{tt},1000);
     fdr_WPM{tt} = pvalue_rank_fdr(WPM_density_pv{tt},1000);
     tmp = [BPM_density_pv{tt} WPM_density_pv{tt}];
     fdr_combined{tt} = pvalue_rank_fdr(tmp,1000);
end

save(sprintf('results_pbody_%s_R%s.mat',model,num2str(R)),'fdr*','BPM*','WPM*')
