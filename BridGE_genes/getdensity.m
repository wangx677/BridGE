function getdensity(ssmFile,R)

% this script compute gene-gene interaction density
% and use ranksum test (-log10 pv) to evaluate interaction based on local background

load BPMind.mat

outputFile = sprintf('genstats_%s_R%s.mat',ssmFile,num2str(R));

% if outputFile doesn't exist
if exist(outputFile, 'file') ~= 2

     load(sprintf('%s_R%s.mat',ssmFile,num2str(R)));
     for tt=1:2
          ssM{tt} = squareform(ssM{tt});
     end

     for tt=1:2
          BPM_density{tt}(1,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM{tt}(x,y),1,[])),BPM.ind1,BPM.ind2,'uniform',0));
          WPM_density{tt}(1,:) = cell2mat(cellfun(@(x,y)mean(reshape(ssM{tt}(x,y),1,[])),WPM.ind,WPM.ind,'uniform',0));

          % use ranksum test to evaluate the signficance of density by comparing the observed GI in BPM/WPM with the pathway background
          bpm_local_1{tt}(1,:) = zeros(size(BPM_density{tt}(1,:)));
          bpm_local_2{tt}(1,:) = zeros(size(BPM_density{tt}(1,:)));

          ind = find(BPM.ind1size>2 & BPM.ind2size>2); % some gene pairs mapped to same SNPs, make sure each gene have at least 2 SNPs are different
          bpm_local_1{tt}(1,ind) = -log10(cell2mat(cellfun(@(x,y)ranksum(reshape(ssM{tt}(x,y),1,[]),reshape(ssM{tt}(x,:),1,[]),'tail','right'),BPM.ind1(ind),BPM.ind2(ind),'uniform',0))); 
          bpm_local_2{tt}(1,ind) = -log10(cell2mat(cellfun(@(x,y)ranksum(reshape(ssM{tt}(x,y),1,[]),reshape(ssM{tt}(y,:),1,[]),'tail','right'),BPM.ind1(ind),BPM.ind2(ind),'uniform',0)));
          bpm_local{tt}(1,:) = min(bpm_local_1{tt}(1,:),bpm_local_2{tt}(1,:));
          wpm_local{tt}(1,:) = -log10(cell2mat(cellfun(@(x)ranksum(reshape(ssM{tt}(x,x),1,[]),reshape(ssM{tt}(x,:),1,[]),'tail','right'),WPM.ind,'uniform',0)));
     end

     save(outputFile,'*density','bpm_local','wpm_local')
end
