function post_LD_filter(genstatsfile)

% load data
load(genstatsfile)
load('BPMind.mat')
load('SNPdataAR.mat')
ssMfile = strsplit(genstatsfile,'genstats_');
ssMfile = ssMfile{2};
load(ssMfile);

ld = readtable('plink.ld','filetype','text');
n = length(SNPdata.rsid);
ldtmp = zeros(n,n);

for i=1:size(ld,1)
     idx1 = find(ismember(SNPdata.rsid,ld.SNP_A{i}));
     idx2 = find(ismember(SNPdata.rsid,ld.SNP_B{i}));
     ldtmp(idx1,idx2) = ld.R2(i);
end

ldtmp = max(ldtmp,ldtmp');


% convert vector to matrix
for tt=1:2
     ssM{tt} = squareform(ssM{tt});
end

% compute SNP-SNP similarity based on interaction
for tt=1:2
     ssMr2{tt} = squareform(1-pdist(ssM{tt}>0.5,'jaccard'));
end

% set GI=0 for redundant SNPs in LD
for tt=1:2
     ssM_ld{tt} = ssM{tt};
     ssMdegree{tt} = sum(ssM{tt}>=0.5);
end

for m=1:length(WPM.ind)
     for tt=1:2
          [tmp idx] = sort(ssMdegree{tt}(WPM.ind{m}),'descend');
          for i=1:length(tmp)-1
               for j=i+1:length(tmp)
                    if ldtmp(WPM.ind{m}(idx(i)),WPM.ind{m}(idx(j)))>0  & ssMr2{tt}(WPM.ind{m}(idx(i)),WPM.ind{m}(idx(j)))>0.75 & SNPdata.chr(WPM.ind{m}(idx(i)))==SNPdata.chr(WPM.ind{m}(idx(j))) && abs(SNPdata.loc(WPM.ind{m}(idx(i)))-SNPdata.loc(WPM.ind{m}(idx(j))))<=300000
                         ssM_ld{tt}(WPM.ind{m}(idx(j)),:) = 0;
                    end
               end
          end
     end
end

for tt=1:2
     ssM_ld{tt} = min(ssM_ld{tt},ssM_ld{tt}');
end

% compute BPM/WPM density and ranksum test 
for tt=1:2
     BPM_density_ld{tt}(1,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM_ld{tt}(x,y),1,[])),BPM.ind1,BPM.ind2,'uniform',0));
     WPM_density_ld{tt}(1,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM_ld{tt}(x,y),1,[])),WPM.ind,WPM.ind,'uniform',0));

     % use ranksum test to evaluate the signficance of density by comparing the observed GI in BPM/WPM with the pathway background
     bpm_local_1{tt}(1,:) = zeros(size(BPM_density_ld{tt}(1,:)));
     bpm_local_2{tt}(1,:) = zeros(size(BPM_density_ld{tt}(1,:)));

     ind = find(BPM.ind1size>2 & BPM.ind2size>2); % some gene pairs mapped to same SNPs, make sure each gene have at least 2 SNPs are different
     bpm_local_1{tt}(1,ind) = -log10(cell2mat(cellfun(@(x,y)ranksum(reshape(ssM_ld{tt}(x,y),1,[]),reshape(ssM_ld{tt}(x,:),1,[]),'tail','right'),BPM.ind1(ind),BPM.ind2(ind),'uniform',0)));
     bpm_local_2{tt}(1,ind) = -log10(cell2mat(cellfun(@(x,y)ranksum(reshape(ssM_ld{tt}(x,y),1,[]),reshape(ssM_ld{tt}(y,:),1,[]),'tail','right'),BPM.ind1(ind),BPM.ind2(ind),'uniform',0)));
     bpm_local_ld{tt}(1,:) = min(bpm_local_1{tt}(1,:),bpm_local_2{tt}(1,:));
     wpm_local_ld{tt}(1,:) = -log10(cell2mat(cellfun(@(x)ranksum(reshape(ssM_ld{tt}(x,x),1,[]),reshape(ssM{tt}(x,:),1,[]),'tail','right'),WPM.ind,'uniform',0)));
end

save(genstatsfile,'BPM_density_ld','WPM_density_ld','bpm_local_ld','wpm_local_ld','ssM_ld','-append')
