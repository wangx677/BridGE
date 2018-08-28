function [bpm_local bpm_local_pv density_bpm density_bpm_expected wpm_local wpm_local_pv density_wpm density_wpm_expected denseidx path_degree path_degree_pv] = genstats(ssmFile,bpmindFile,binaryNetwork,snpPerms,minPath,netDensity)

% FUNCTION [bpm_local density_bpm density_bpm_expected wpm_local density_wpm density_wpm_expected denseidx path_degree] = genstats(ssmFile,bpmindFile,minPath,binaryNetwork,netDensity) 
%
% GENSTATS calculates general statistics for BPM (chi-square global,chi-square local, density), WPM (chi-square global,chi-square local, density) and 
%    PATH (ranksum)
%
% INPUTS:
%   ssmFile - SNP-SNP interaction file name
%   bpmindFile - a mat-file that includes SNP indexes of pathways for all possbile BPMs, WPMs
%   binaryNetwork - 1: binary network; 0: weighted network
%   snpPerms - number of SNP permutation
%   netDensity - network binarization density based on considering both protecitve/risk interactions
% 
% OUTPUTS:
%   An output file named genstats_<ssmFile>.mat or genstats_<ssmFile>_density<netDensity>.mat
%   

if (nargin > 6)
    error('Incorrect number of input arguments!');
end

load(bpmindFile)

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
   [bpm_local{tt} bpm_local_pv{tt} density_bpm{tt} density_bpm_expected{tt} wpm_local{tt} wpm_local_pv{tt} density_wpm{tt} density_wpm_expected{tt} denseidx{tt} path_degree{tt} path_degree_pv{tt}] = rungenstats(full(ssM{tt}),BPM.ind1,BPM.ind2,BPM.ind1size,BPM.ind2size,WPM.ind,bpmsize,wpmsize,pathsize,binaryNetwork,snpPerms,minPath);
end

outputFile = sprintf('genstats_%s.mat',ssmFile);

if binaryNetwork==1
     save(outputFile,'bpm_local*','wpm_local*','path_degree*','denseidx','density_*','netDensity','-v7.3')
else
     save(outputFile,'bpm_local*','wpm_local*','path_degree*','denseidx','density_*','-v7.3')
end

end

function [bpm_local bpm_local_pv density_bpm density_bpm_expected wpm_local wpm_local_pv density_wpm density_wpm_expected denseidx path_degree path_degree_pv] = rungenstats(MM,BPMind1,BPMind2,BPMind1size,BPMind2size,WPMind,bpmsize,wpmsize,pathsize,binaryNetwork,snpPerms,minPath) 

[p q] = size(MM);
N = p;

% BPM
% For weighted network, first filter BPMs using binarized network since the network is very sparse
% only BPMs with chi2_local p<=0.05 will be further evaluated for bpm_local with ranksum test with weighted network.
% For binary network, bpm_local is based on chi2 test.
tic
% chi2 BPM local
if  binaryNetwork==1 & sum(sum(MM>1))==0 % confirm MM is a binary network
     % get interaction counts in BPM for binary network 
     bpmgi = zeros(1,length(BPMind1));
     for i=1:length(BPMind1);bpmgi(i) = sum2(MM(BPMind1{i},BPMind2{i}));end
     bpmnotgi = bpmsize-bpmgi;
     density_bpm = bpmgi./bpmsize;
     density_bpm(isnan(density_bpm)==1) = 0;

     % get interaction counts in pathway subnetwork for binary network
     sumMM = sum(MM,2);
     path1bggi = zeros(1,length(BPMind1));
     path2bggi = zeros(1,length(BPMind2));
     for i=1:length(BPMind1);path1bggi(i) = sum(sumMM(BPMind1{i}));end
     for i=1:length(BPMind2);path2bggi(i) = sum(sumMM(BPMind2{i}));end
     
     path1bggi = path1bggi-bpmgi;% pathway1 background interaction counts
     path2bggi = path2bggi-bpmgi;% pathway2 background interaction counts

elseif binaryNetwork==0 & sum(sum(abs(MM)>1))>0 % confirm MM is a weighted network    
     % for filtering purpose in weighted network
     % 0.2 is the cutoff for hygeSSI network to avoid super low hygeSSI scores
     MMtmp = MM>=0.2;
     bpmgi = zeros(1,length(BPMind1));
     for i=1:length(BPMind1);bpmgi(i) = sum2(MMtmp(BPMind1{i},BPMind2{i}));end
     bpmnotgi = bpmsize-bpmgi;

     % get interaction counts in pathway subnetwork for weighted network
     sumMM = sum(MMtmp,2);
     path1bggi = zeros(1,length(BPMind1));
     path2bggi = zeros(1,length(BPMind2));
     for i=1:length(BPMind1);path1bggi(i) = sum(sumMM(BPMind1{i}));end
     for i=1:length(BPMind2);path2bggi(i) = sum(sumMM(BPMind2{i}));end

     path1bggi = path1bggi-bpmgi; % pathway1 background interaction counts
     path2bggi = path2bggi-bpmgi; % pathway2 background interaction counts

else
     warning('SNP-SNP interaction network is inconsistent with network type.')
     return  
end


path1bgsize = BPMind1size*p; % pathway1 background size
path2bgsize = BPMind2size*p; % pathway2 background size
path1bgnotgi = (path1bgsize-bpmsize)-path1bggi; % pathway1 background not interaction counts
path2bgnotgi = (path2bgsize-bpmsize)-path2bggi; % pathway2 background not interaction counts

% use chi2-square test to check BPM enrichment significance comparing to two local pathways
chi2_bpm_local_1 = -log10(1-chi2cdf(chi2test([bpmgi' path1bggi' bpmnotgi' path1bgnotgi']),1));
chi2_bpm_local_1(bpmgi./(bpmgi+bpmnotgi) < path1bggi./(path1bggi+path1bgnotgi)) = -chi2_bpm_local_1(bpmgi./(bpmgi+bpmnotgi) < path1bggi./(path1bggi+path1bgnotgi));

chi2_bpm_local_2 = -log10(1-chi2cdf(chi2test([bpmgi' path2bggi' bpmnotgi' path2bgnotgi']),1));
chi2_bpm_local_2(bpmgi./(bpmgi+bpmnotgi) < path2bggi./(path2bggi+path2bgnotgi)) = -chi2_bpm_local_2(bpmgi./(bpmgi+bpmnotgi) < path2bggi./(path2bggi+path2bgnotgi));

density_bpm_local_1 = ((bpmgi+path1bggi)./(path1bgnotgi+path1bggi+bpmsize))';
density_bpm_local_2 = ((bpmgi+path2bggi)./(path2bgnotgi+path2bggi+bpmsize))';

clear bpmgi path1bggi bpmnotgi path1bgnotgi path2bggi path2bgnotgi

toc1 = toc
sprintf('BPM chi2 walltime: %f',toc1)

% We will pick the denser pathway (less significant chi2 local score) for rank sum test.
% Since we only care about enrichment of BPM interactions with respect to local pathways,
% we will do SNP permutation on one pathway (the denser one) dimention.

% Each BPM has BPMind1 and BPMind2 that are two sets of pathway SNP index,
% and we need to find out which pathway to be used as SNP permutation background.
% To do so, we swtich BPMind1 and BPMind2 based on their densities or chi2 local scores,
% so that we only need to shuffle one side of the SNP-SNP interaction matrix.

% denseidx points to the set of denser pathways and the corresponding tests are less (min) signficant
% chi2_bpm_local from the denser pathways will be used as final score

tic
% need to consider different scenarios
density_bpm_local_1(isnan(density_bpm_local_1)) = 0;
density_bpm_local_2(isnan(density_bpm_local_2)) = 0;

denseidx = zeros(length(density_bpm_local_1),1);

% scenarios 1: chi2_local tests agree with density
maxind1 = find(chi2_bpm_local_1<chi2_bpm_local_2 & density_bpm_local_1>density_bpm_local_2);
maxind2 = find(chi2_bpm_local_1>chi2_bpm_local_2 & density_bpm_local_1<density_bpm_local_2);

denseidx(maxind1) = 1;
denseidx(maxind2) = 2;

ind = find(denseidx == 0);

% scenarios 2: chi2_local tests are equal, but densities are not
maxind1 = find(chi2_bpm_local_1(ind)==chi2_bpm_local_2(ind) & density_bpm_local_1(ind)>density_bpm_local_2(ind));
maxind2 = find(chi2_bpm_local_1(ind)==chi2_bpm_local_2(ind) & density_bpm_local_1(ind)<density_bpm_local_2(ind));

denseidx(ind(maxind1)) = 1;
denseidx(ind(maxind2)) = 2;

ind = find(denseidx == 0);

% scenarios 3: densities are equal, but chi2_local tests are not
maxind1 = find(chi2_bpm_local_1(ind)<chi2_bpm_local_2(ind) & density_bpm_local_1(ind)==density_bpm_local_2(ind));
maxind2 = find(chi2_bpm_local_1(ind)>chi2_bpm_local_2(ind) & density_bpm_local_1(ind)==density_bpm_local_2(ind));

denseidx(ind(maxind1)) = 1;
denseidx(ind(maxind2)) = 2;

ind = find(denseidx == 0);

% scenarios 4: chi2_local and densities don't agree with each other and are not equal, use chi2_local to decide densier pathway
maxind1 = find(chi2_bpm_local_1(ind)<=chi2_bpm_local_2(ind));
maxind2 = find(chi2_bpm_local_1(ind)>chi2_bpm_local_2(ind));

denseidx(ind(maxind1)) = 1;
denseidx(ind(maxind2)) = 2;

clear maxind1 maxind2

% finalize chi2_local scores
chi2_bpm_local = zeros(1,length(chi2_bpm_local_1));
chi2_bpm_local(find(denseidx==1)) = chi2_bpm_local_1(find(denseidx==1));;
chi2_bpm_local(find(denseidx==2)) = chi2_bpm_local_2(find(denseidx==2));

clear chi2_bpm_local_1 chi2_bpm_local_2 

% only evaluate  BPMs with chi2_bpm_local p<=0.05 and 
% the number of SNPs in each size of BPM is greater than minPath for further analysis
ind2keep_bpm = find(chi2_bpm_local>=-log10(0.1) & BPMind1size>=minPath & BPMind2size>=minPath);

% swap BPMind1 and BPMind2
% so that BPMind1new is for pathway with denser subnetwork
% BPMind2new is for pathway with sparser network
% so for network MM(x,y),the SNP shuffling will be on y dimention 
BPMind1new = cell(1,length(BPMind1));
BPMind2new = cell(1,length(BPMind2));

BPMind1new(find(denseidx==1)) = BPMind1(find(denseidx==1));
BPMind1new(find(denseidx==2)) = BPMind2(find(denseidx==2));

BPMind2new(find(denseidx==1)) = BPMind2(find(denseidx==1));
BPMind2new(find(denseidx==2)) = BPMind1(find(denseidx==2));

toc2 = toc

sprintf('BPM denser pathway walltime: %f',toc2)

tic
% WPM local
if binaryNetwork==1 & sum(sum(MM>1))==0 % confirm MM is a binary network
     % get interaction counts in WPM for binary network 
     wpmgi = zeros(1,length(WPMind));
     for i=1:length(WPMind);wpmgi(i) = sum2(MM(WPMind{i},WPMind{i}));end
     wpmnotgi = wpmsize-wpmgi;
     density_wpm = wpmgi./wpmsize;
     density_wpm(isnan(density_wpm)==1) = 0;

     % get interaction counts in pathway subnetwork for binary network
     sumMM = sum(MM,2);
     pathbggi = zeros(1,length(WPMind));
     for i=1:length(WPMind);pathbggi(i) = sum(sumMM(WPMind{i}));end
     
     pathbggi = pathbggi-wpmgi;% pathway background interaction counts

elseif binaryNetwork==0 & sum(sum(abs(MM)>1))>0 % confirm MM is a weighted network
     % for filtering purpose in weighted network
     % 0.2 is the cutoff for hygeSSI network to avoid super low hygeSSI scores
     MMtmp = MM>=0.2;
     wpmgi = zeros(1,length(WPMind));
     for i=1:length(WPMind);wpmgi(i) = sum2(MMtmp(WPMind{i},WPMind{i}));end
     wpmnotgi = wpmsize-wpmgi;

     % get interaction counts in pathway subnetwork for binary network
     sumMM = sum(MMtmp,2);
     pathbggi = zeros(1,length(WPMind));
     for i=1:length(WPMind);pathbggi(i) = sum(sumMM(WPMind{i}));end

     pathbggi = pathbggi-wpmgi;% pathway background interaction counts

end

pathbgsize = pathsize*p; % pathway background size
pathbgnotgi = (pathbgsize-wpmsize)-pathbggi; % pathway background not interaction counts

chi2_wpm_local = -log10(1-chi2cdf(chi2test([wpmgi' pathbggi' wpmnotgi' pathbgnotgi']),1));
chi2_wpm_local(wpmgi./(wpmgi+wpmnotgi) < pathbggi./(pathbggi+pathbgnotgi)) = -chi2_wpm_local(wpmgi./(wpmgi+wpmnotgi) < pathbggi./(pathbggi+pathbgnotgi));

ind2keep_wpm = find(chi2_wpm_local>=-log10(0.1));

clear wpmgi pathbggi wpmnotgi pathbgnotgi

toc3 = toc
sprintf('WPM chi2 walltime: %f',toc3)


% BPM WPM ranksum
% compute BPM and WPM scores based on ranksum test if it is a weighted network
tic
if binaryNetwork==0  
     density_bpm = zeros(size(BPMind1new));
     bpm_local = zeros(size(BPMind1new));
     bpmsum = zeros(size(BPMind1new));

     ind1 = BPMind1new(ind2keep_bpm);
     ind2 = BPMind2new(ind2keep_bpm);
     
     bpmsum_tmp = zeros(1,length(ind1));
     for i=1:length(ind1); bpmsum_tmp(i) = sum2(MM(ind1{i},ind2{i}));end
     density_bpm_tmp = bpmsum_tmp./bpmsize(ind2keep_bpm);

     % modified ranksum test for efficiency
     
     bpm_local_tmp = zeros(1,length(ind2));
     for i=1:length(ind1)
          MMtmp = MM(ind1{i},:);
          bpm_local_tmp(i) = rs(MMtmp,ind2{i});
     end

     density_bpm(ind2keep_bpm) = density_bpm_tmp;
     bpm_local(ind2keep_bpm) = -log10(bpm_local_tmp);
     bpmsum(ind2keep_bpm) = bpmsum_tmp;
     clear density_bpm_tmp bpm_local_tmp ind1 ind2

     % WPM
     density_wpm = zeros(size(WPMind));
     wpm_local = zeros(size(WPMind));
     wpmsum = zeros(size(WPMind));

     ind1 = WPMind(ind2keep_wpm);
     wpmsum_tmp = cellfun(@(x)sum2(MM(x,x)),ind1);
     density_wpm_tmp = wpmsum_tmp./wpmsize(ind2keep_wpm);
     
     wpm_local_tmp = zeros(1,length(ind1));
     for i=1:length(ind1)
          MMtmp = MM(ind1{i},:);
          wpm_local_tmp(i) = rs(MMtmp,ind1{i});
     end
     
     density_wpm(ind2keep_wpm) = density_wpm_tmp;
     wpm_local(ind2keep_wpm) = -log10(wpm_local_tmp);
     wpmsum(ind2keep_wpm) = wpmsum_tmp;

     ind2keep_bpm = find(bpm_local>=-log10(0.05));
     ind2keep_wpm = find(wpm_local>=-log10(0.05));
     clear wpmsum_tmp density_wpm_tmp wpm_local_tmp ind1  
else
     bpm_local = chi2_bpm_local;
     wpm_local = chi2_wpm_local;
     
     % BPM
     density_bpm = zeros(size(BPMind1new));
     bpmsum = zeros(size(BPMind1new));

     ind1 = BPMind1new(ind2keep_bpm);
     ind2 = BPMind2new(ind2keep_bpm);

     bpmsum_tmp = zeros(1,length(ind1));
     for i=1:length(ind1); bpmsum_tmp(i) = sum2(MM(ind1{i},ind2{i}));end
     density_bpm_tmp = bpmsum_tmp./bpmsize(ind2keep_bpm);

     density_bpm(ind2keep_bpm) = density_bpm_tmp;
     bpmsum(ind2keep_bpm) = bpmsum_tmp;
     clear density_bpm_tmp ind1 ind2

     % WPM
     density_wpm = zeros(size(WPMind));
     wpmsum = zeros(size(WPMind));

     ind1 = WPMind(ind2keep_wpm);
     wpmsum_tmp = cellfun(@(x)sum2(MM(x,x)),ind1);
     density_wpm_tmp = wpmsum_tmp./wpmsize(ind2keep_wpm);

     density_wpm(ind2keep_wpm) = density_wpm_tmp;
     wpmsum(ind2keep_wpm) = wpmsum_tmp;

     clear wpmsum_tmp density_wpm_tmp ind1
end
toc4 = toc

sprintf('BPM WPM ranksum walltime: %f',toc4)

% compute expected BPM and WPM density
tic
density_bpm_expected = zeros(size(BPMind1));

ind1 = BPMind1new(ind2keep_bpm);
sumMM = sum(MM);

for i=1:length(BPMind1new);density_bpm_expected(i) = sum(sumMM(BPMind1new{i}))/(N*length(BPMind1new{i}));end

density_wpm_expected = zeros(size(WPMind));

for i=1:length(WPMind);density_wpm_expected(i) = sum(sumMM(WPMind{i}))/(N*length(WPMind{i}));end

toc5 = toc
sprintf('BPM WPM expected density walltime: %f',toc5)

% PATH
path_degree = ones(1,length(WPMind));
for i=1:length(WPMind);path_degree(i) = ranksum(sumMM(WPMind{i}),sumMM(setdiff(1:N,WPMind{i})),'tail','right');end
path_degree = -log10(path_degree);
ind2keep_path = find(path_degree>=-log10(0.1));

% launch snpPerms permutation
% need to consider ind2keep_* is empty
tic
count_bpm = zeros(size(ind2keep_bpm));
count_wpm = zeros(size(ind2keep_wpm));
count_path = zeros(size(ind2keep_path));

bpmind1 = BPMind1new(ind2keep_bpm);
bpmind2 = BPMind2new(ind2keep_bpm);
wpmind = WPMind(ind2keep_wpm);
pathind = WPMind(ind2keep_path);

bpm_local_pv = ones(size(BPMind1));
wpm_local_pv = ones(size(WPMind));
path_degree_pv = ones(size(WPMind));

for j=1:snpPerms
     rand('seed',j)
     MMtmp = MM(:,randperm(N));
     sumMMtmp = sum(MMtmp);
     
     bpmsum_tmp = zeros(1,length(bpmind1));
     for i=1:length(bpmind1);bpmsum_tmp(i) = sum2(MMtmp(bpmind1{i},bpmind2{i}));end 
     count_bpm = reshape(count_bpm,1,length(bpmind1))+reshape((bpmsum_tmp >= bpmsum(ind2keep_bpm)),1,length(bpmind1));

     wpmsum_tmp = zeros(1,length(wpmind));
     for i=1:length(wpmind);wpmsum_tmp(i)=sum2(MMtmp(wpmind{i},wpmind{i}));end
     count_wpm = reshape(count_wpm,1,length(wpmind))+reshape((wpmsum_tmp >= wpmsum(ind2keep_wpm)),1,length(wpmind));

     path_degree_tmp = zeros(1,length(pathind));
     for i=1:length(pathind);path_degree_tmp(i) = ranksum(sumMMtmp(pathind{i}),sumMMtmp(setdiff(1:N,pathind{i})));end
     count_path = reshape(count_path,1,length(pathind))+reshape((path_degree_tmp>=path_degree(ind2keep_path)),1,length(pathind));
     
     clear bpmsum_tmp wpmsum_tmp path_degree_tmp 
end 

bpm_local_pv(ind2keep_bpm) = (count_bpm+1)/snpPerms;
wpm_local_pv(ind2keep_wpm) = (count_wpm+1)/snpPerms;
path_degree_pv(ind2keep_path) = (count_path+1)/snpPerms;
toc6 = toc
disp(sprintf('SNP permutation costs %s',toc6))

end
