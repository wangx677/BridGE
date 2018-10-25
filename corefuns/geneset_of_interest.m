function geneset_of_interest(genes,ssmFile,snpgenemapping,outputname,binaryNetwork,snpPerms,minPath,netDensity)

% this script is used to identify pathways that are interacting with the given geneset

load(genes)
load(snpgenemapping)

% get snp index that can be mapped to genes
ind = find(ismember(snp2gene.genelist,genes));
ind = find(sum(snp2gene.sgmatrix(:,ind),2)~=0);

load SNPdataAR.mat
targt_snps = SNPdata.rsid(ind);

% get BPM infomation
load BPMind.mat
clear BPM
for i=1:length(WPM.ind)
     ind1{i} = ind(find(ismember(ind,WPM.ind{i})==0));
     ind2{i} = WPM.ind{i}(find(ismember(WPM.ind{i},ind)==0));
end

BPM.ind1 = [ind ind1];
BPM.ind2 = [ind ind2];
BPM.ind1size = cellfun(@(x)length(x),BPM.ind1);
BPM.ind2size = cellfun(@(x)length(x),BPM.ind2);
BPM.size = BPM.ind1size.*BPM.ind2size;
save(sprintf('BPMind_%s.mat',outputname),'BPM')

% calculate stats
load(sprintf('%s.mat',ssmFile))
[p q] = size(ssM{1});
if min(p,q) == 1;
     for tt=1:2
          ssM{tt} = squareform(ssM{tt});
     end
end

[p q] = size(ssM{1});
if binaryNetwork==1 & exist('netDensity')==1
     tmp = max(ssM{1},ssM{2});
     tmpcut = quantile(tmp(:),1-netDensity);
     for tt=1:2
          ssM{tt} = ssM{tt}>=tmpcut;
     end

     clear tmp tmpcut
elseif binaryNetwork==1 & exist('netDensity')~=1
     for tt=1:2
          ssM{tt} = ssM{tt}>0;
     end
     netDensity = nnz(max(ssM{1},ssM{2}))/(p*q)
end

% analyze protective and risk effect individually
bpmsize = BPM.size;
wpmsize = WPM.size;
pathsize = WPM.indsize;
for tt=1:2
   [bpm_local{tt} bpm_local_pv{tt} density_bpm{tt} density_bpm_expected{tt} denseidx{tt}] = rungenstats(full(ssM{tt}),BPM.ind1,BPM.ind2,BPM.ind1size,BPM.ind2size,BPM.size,binaryNetwork,snpPerms,minPath);
end

outputFile = sprintf('genstats_%s_%s.mat',outputname,ssmFile);

if binaryNetwork==1
     save(outputFile,'bpm_local*','denseidx','density_*','netDensity','genes','-v7.3')
else
     save(outputFile,'bpm_local*','denseidx','density_*','genes','-v7.3')
end

