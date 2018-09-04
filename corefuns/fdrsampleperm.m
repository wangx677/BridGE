function fdrsampleperm(model,interaction,pcut,minPath,N)

% model='combined';
% interaction
% pcut:
% N:

load BPMind.mat
BPMsize = repmat(reshape(BPM.size,1,length(BPM.size)),1,2);
WPMsize = repmat(reshape(WPM.size,1,length(WPM.size)),1,2);

BPMpath1size = repmat(reshape(BPM.ind1size,1,length(BPM.ind1size)),1,2);
BPMpath2size = repmat(reshape(BPM.ind2size,1,length(BPM.ind2size)),1,2);

clear BPM WPM

for i=0:N
     if strcmp(interaction,'hygeSSI')
          load(sprintf('genstats_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R%s.mat',model,num2str(i)))
     elseif strcmp(interaction,'regSSI')
          load(sprintf('genstats_ssM_regSSI_%s_R%s.mat',model,num2str(i)))
     elseif strcmp(interaction,'lrSSI')
          load(sprintf('ssM_lr_cassi_pv0.05_R%s.mat',num2str(i)))
     end
     bpm(i+1,:) = [bpm_local{1} bpm_local{2}];
     wpm(i+1,:) = [wpm_local{1} wpm_local{2}];
     path(i+1,:) = [path_degree{1} path_degree{2}];
     
     bpm_pv(i+1,:) = [bpm_local_pv{1} bpm_local_pv{2}];
     wpm_pv(i+1,:) = [wpm_local_pv{1} wpm_local_pv{2}];
     path_pv(i+1,:) = [path_degree_pv{1} path_degree_pv{2}];
end

clear bpm_local wpm_local bpm_local_pv wpm_local_pv path_degree path_degree_pv density* maxidx minidx 

IND = find(BPMpath1size<minPath | BPMpath2size<minPath);
bpm(:,IND) = 0;
bpm_pv(:,IND) = 1;

%% compute FDR for BPMs
ind = find(bpm_pv(1,:)<=pcut);

if isempty(ind)~=1

for i=1:length(ind)
     % based on local ranksum test p-value
     m1(i) = nnz(bpm_pv(1,:)<=bpm_pv(1,ind(i)));
     m2(i) = nnz(bpm_pv(2:end,:)<=bpm_pv(1,ind(i)))/N;

     % combine local ranksum test and its p-value
     n1(i) = nnz(bpm_pv(1,:)<=bpm_pv(1,ind(i)) & bpm(1,:)>=round(bpm(1,ind(i))));
     n2(i) = nnz(bpm_pv(2:end,:)<=bpm_pv(1,ind(i)) & bpm(2:end,:)>= round(bpm(1,ind(i))))/N;
end

fdr1 = m2./m1;
fdr2 = n2./n1;

clear m1 m2 n1 n2
 
fdrBPM1 = ones(size(bpm(1,:)));
fdrBPM2 = ones(size(bpm(1,:)));

% sort fdr based on bpm_pv scores 
testM = [fdr1' bpm_pv(1,ind)'];
[testM idx] = sortrows(testM,[-2]);

for i=1:size(testM,1)
     testM(i,1) = min(testM(1:i,1));
end

% assign FDR to BPMs
fdrBPM1(ind(idx)) = testM(:,1);

% sort fdr based on bpm_pv and bpm_degree_local
testM = [fdr2' bpm_pv(1,ind)' round(bpm(1,ind)')];

[testM idx] = sortrows(testM,[-2 3]);

for i=1:size(testM,1)
     testM(i,1) = min(testM(1:i,1));
end

% assign FDR to BPMs
fdrBPM2(ind(idx)) = testM(:,1);
end


%% compute FDR for WPMs
ind = find(wpm_pv(1,:)<=0.05);

if isempty(ind)~=1

for i=1:length(ind)
     % based on local ranksum test p-value
     m1(i) = nnz(wpm_pv(1,:)<=wpm_pv(1,ind(i)));
     m2(i) = nnz(wpm_pv(2:end,:)<=wpm_pv(1,ind(i)))/N;

     % combine local ranksum test and its p-value
     n1(i) = nnz(wpm_pv(1,:)<=wpm_pv(1,ind(i)) & wpm(1,:)>=wpm(1,ind(i)));
     n2(i) = nnz(wpm_pv(2:end,:)<=wpm_pv(1,ind(i)) & wpm(2:end,:)>= wpm(1,ind(i)))/N;
end

fdr1 = m2./m1;
fdr2 = n2./n1;

clear m1 m2 n1 n2

fdrWPM1 = ones(size(wpm(1,:)));
fdrWPM2 = ones(size(wpm(1,:)));

% sort fdr based on wpm_pv scores
testM = [fdr1' wpm_pv(1,ind)'];
[testM idx] = sortrows(testM,[-2]);

for i=1:size(testM,1)
     testM(i,1) = min(testM(1:i,1));
end

% assign FDR to WPMs
fdrWPM1(ind(idx)) = testM(:,1);

% sort fdr based on wpm_pv and wpm_degree_local
testM = [fdr2' wpm_pv(1,ind)' round(wpm(1,ind)')];

[testM idx] = sortrows(testM,[-2 3]);

for i=1:size(testM,1)
     testM(i,1) = min(testM(1:i,1));
end

% assign FDR to WPMs
fdrWPM2(ind(idx)) = testM(:,1);
end

%% compute FDR for PATH

ind = find(path_pv(1,:)<=0.05);
if isempty(ind)~=1
 
for i=1:length(ind)
     % based on geometric mean of local ranksum test
     m1(i) = nnz(path_pv(1,:)<=path_pv(1,ind(i)));
     m2(i) = nnz(path_pv(2:end,:)<=path_pv(1,ind(i)))/N;
    
     n1(i) = nnz(path_pv(1,:)<=path_pv(1,ind(i)) & path(1,:)>=path(1,ind(i)));
     n2(i) = nnz(path_pv(2:end,:)<=path_pv(1,ind(i)) & path(2:end,:)>=path(1,ind(i)))/N; 
end

fdr1 = m2./m1;
fdr2 = n2./n1;

clear m1 m2 n1 n2
 
fdrPATH1 = ones(size(path(1,:)));
fdrPATH2 = ones(size(path(1,:)));

% sort fdr based on bpm scores 
testM = [fdr1' path_pv(1,ind)'];
[testM idx] = sortrows(testM,[-2]);

for i=1:size(testM,1)
     testM(i,1) = min(testM(1:i,1));
end

% assign FDR to BPMs
fdrPATH1(ind(idx)) = testM(:,1);

% sort fdr based on bpm scores
testM = [fdr2' path_pv(1,ind)' path(1,ind)'];
[testM idx] = sortrows(testM,[-2 3]);

for i=1:size(testM,1)
     testM(i,1) = min(testM(1:i,1));
end

clear fdr1 fdr2

% assign FDR to BPMs
fdrPATH2(ind(idx)) = testM(:,1);
end
bpm_ranksum = bpm(1,:);
wpm_ranksum = wpm(1,:);
path_ranksum = path(1,:);

bpm_pv = bpm_pv(1,:);
wpm_pv = wpm_pv(1,:);
path_pv = path_pv(1,:);
if ismember(interaction,'hygeSSI')==1
     save(sprintf('results_ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',model),'fdr*','*_ranksum','bpm_pv','wpm_pv','path_pv')
elseif ismember(interaction,'regSSI')==1
     save(sprintf('results_ssM_regSSI_%s_R0.mat',model),'fdr*','*_ranksum','bpm_pv','wpm_pv','path_pv')
end
