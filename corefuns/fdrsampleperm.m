function fdrsampleperm(ssmFile,BPMindFile,pcut,minPath,N,genesetname)

% model='combined';
% interaction
% pcut:
% N:

load(BPMindFile)
BPMsize = repmat(reshape(BPM.size,1,length(BPM.size)),1,2);

% in case WPM is not evaluated
if exist('WPM','var')
     WPMsize = repmat(reshape(WPM.size,1,length(WPM.size)),1,2);
     clear WPM
end

BPMpath1size = repmat(reshape(BPM.ind1size,1,length(BPM.ind1size)),1,2);
BPMpath2size = repmat(reshape(BPM.ind2size,1,length(BPM.ind2size)),1,2);

clear BPM

for i=0:N
     if exist('genesetname','var')==0
          genstatsfile = sprintf('genstats_%s_R%s.mat',ssmFile,num2str(i));
     else
          genstatsfile = sprintf('genstats_%s_%s_R%s.mat',genesetname,ssmFile,num2str(i));
     end
     load(genstatsfile)
     bpm(i+1,:) = [bpm_local{1} bpm_local{2}];
     bpm_pv(i+1,:) = [bpm_local_pv{1} bpm_local_pv{2}];

     bpm_protective(i+1,:) = bpm_local{1};
     bpm_risk(i+1,:) =  bpm_local{2};

     bpm_pv_protective(i+1,:) = bpm_local_pv{1};
     bpm_pv_risk(i+1,:) =  bpm_local_pv{2};

     if exist('WPMsize','var')
          wpm(i+1,:) = [wpm_local{1} wpm_local{2}];
          wpm_protective(i+1,:) = wpm_local{1};
          wpm_risk(i+1,:) =  wpm_local{2};

          path(i+1,:) = [path_degree{1} path_degree{2}];
          path_protective(i+1,:) = path_degree{1};
          path_risk(i+1,:) =  path_degree{2};

          wpm_pv(i+1,:) = [wpm_local_pv{1} wpm_local_pv{2}];
          wpm_pv_protective(i+1,:) = wpm_local_pv{1};
          wpm_pv_risk(i+1,:) =  wpm_local_pv{2};

          path_pv(i+1,:) = [path_degree_pv{1} path_degree_pv{2}];
          path_pv_protective(i+1,:) = path_degree_pv{1};
          path_pv_risk(i+1,:) =  path_degree_pv{2};

          clear wpm_local wpm_local_pv path_degree path_degree_pv
     end

end

clear bpm_local bpm_local_pv density* maxidx minidx 

IND = find(BPMpath1size<minPath | BPMpath2size<minPath);
bpm(:,IND) = 0;
bpm_pv(:,IND) = 1;

%% compute FDR for BPMs (protective and risk)
ind = find(bpm_pv(1,:)<=pcut);

if isempty(ind)~=1

     for i=1:length(ind)
          % combine local ranksum test and its p-value
          n1(i) = nnz(bpm_pv(1,:)<=bpm_pv(1,ind(i)) & bpm(1,:)>=bpm(1,ind(i)));
          n2(i) = nnz(bpm_pv(2:end,:)<=bpm_pv(1,ind(i)) & bpm(2:end,:)>= bpm(1,ind(i)))/N;
     end

     fdr = n2./n1;

     clear n1 n2
 
     fdrBPM = ones(size(bpm(1,:)));

     % adjust FDR
     % sometimes more significant BPMs can have higher FDRs
     % we make FDR adjustment to such BPMs
     testM = [fdr' bpm_pv(1,ind)' bpm(1,ind)'];

     [testM idx] = sortrows(testM,-1);

     for i=1:size(testM,1)
          indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
          fdrnew(i) = min(testM(indtmp,1));
     end

     % assign FDR to BPMs
     fdrBPM(ind(idx)) = fdrnew;

     clear fdr fdrnew 
end

%% compute FDR for BPMs (protective only)
ind = find(bpm_pv_protective(1,:)<=pcut);

if isempty(ind)~=1

     for i=1:length(ind)
          % combine local ranksum test and its p-value
          n1(i) = nnz(bpm_pv_protective(1,:)<=bpm_pv_protective(1,ind(i)) & bpm_protective(1,:)>=bpm_protective(1,ind(i)));
          n2(i) = nnz(bpm_pv_protective(2:end,:)<=bpm_pv_protective(1,ind(i)) & bpm_protective(2:end,:)>= bpm_protective(1,ind(i)))/N;
     end

     fdr = n2./n1;

     clear n1 n2

     fdrBPM_protective = ones(size(bpm_protective(1,:)));

     % sort fdr based on bpm_pv_protective and bpm_degree_local
     testM = [fdr' bpm_pv_protective(1,ind)' bpm_protective(1,ind)'];

     [testM idx] = sortrows(testM,-1);

     for i=1:size(testM,1)
          indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
          fdrnew(i) = min(testM(indtmp,1));
     end

     % assign FDR to BPMs
     fdrBPM_protective(ind(idx)) = fdrnew;

     clear fdr fdrnew
else
     fdrBPM_protective = ones(size(bpm_protective(1,:)));
end

%% compute FDR for BPMs (risk only)
ind = find(bpm_pv_risk(1,:)<=pcut);

if isempty(ind)~=1

     for i=1:length(ind)
          % combine local ranksum test and its p-value
          n1(i) = nnz(bpm_pv_risk(1,:)<=bpm_pv_risk(1,ind(i)) & bpm_risk(1,:)>=bpm_risk(1,ind(i)));
          n2(i) = nnz(bpm_pv_risk(2:end,:)<=bpm_pv_risk(1,ind(i)) & bpm_risk(2:end,:)>= bpm_risk(1,ind(i)))/N;
     end

     fdr = n2./n1;

     clear n1 n2

     fdrBPM_risk = ones(size(bpm_risk(1,:)));

     % sort fdr based on bpm_pv_risk and bpm_degree_local
     testM = [fdr' bpm_pv_risk(1,ind)' bpm_risk(1,ind)'];

     [testM idx] = sortrows(testM,-1);

     for i=1:size(testM,1)
          indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
          fdrnew(i) = min(testM(indtmp,1));
     end

     % assign FDR to BPMs
     fdrBPM_risk(ind(idx)) = fdrnew;

     clear fdr fdrnew
else
     fdrBPM_risk = ones(size(bpm_risk(1,:)));
end


if exist('WPMsize','var')
     %% compute FDR for WPMs (protective and risk)
     ind = find(wpm_pv(1,:)<=0.05);

     if isempty(ind)~=1

          for i=1:length(ind)
               % combine local ranksum test and its p-value
               n1(i) = nnz(wpm_pv(1,:)<=wpm_pv(1,ind(i)) & wpm(1,:)>=wpm(1,ind(i)));
               n2(i) = nnz(wpm_pv(2:end,:)<=wpm_pv(1,ind(i)) & wpm(2:end,:)>= wpm(1,ind(i)))/N;
          end

          fdr = n2./n1;

          clear n1 n2

          fdrWPM = ones(size(wpm(1,:)));

          % sort fdr based on wpm_pv and wpm_degree_local
          testM = [fdr' wpm_pv(1,ind)' round(wpm(1,ind)')];

          [testM idx] = sortrows(testM,-1);

          for i=1:size(testM,1)
               indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
               fdrnew(i) = min(testM(indtmp,1));
          end

          % assign FDR to WPMs
          fdrWPM(ind(idx)) = fdrnew;
          clear fdr fdrnew
     else
          fdrWPM = ones(size(wpm(1,:)));
     end

     %% compute FDR for WPMs (protective only)
     ind = find(wpm_pv_protective(1,:)<=0.05);

     if isempty(ind)~=1

          for i=1:length(ind)
               % combine local ranksum test and its p-value
               n1(i) = nnz(wpm_pv_protective(1,:)<=wpm_pv_protective(1,ind(i)) & wpm_protective(1,:)>=wpm_protective(1,ind(i)));
               n2(i) = nnz(wpm_pv_protective(2:end,:)<=wpm_pv_protective(1,ind(i)) & wpm_protective(2:end,:)>= wpm_protective(1,ind(i)))/N;
          end

          fdr = n2./n1;

          clear n1 n2

          fdrWPM_protective = ones(size(wpm_protective(1,:)));

          % sort fdr based on wpm_pv_protective and wpm_degree_local
          testM = [fdr' wpm_pv_protective(1,ind)' round(wpm_protective(1,ind)')];

          [testM idx] = sortrows(testM,-1);

           for i=1:size(testM,1)
               indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
               fdrnew(i) = min(testM(indtmp,1));
          end

          % assign FDR to WPMs
          fdrWPM_protective(ind(idx)) = fdrnew;
          clear fdr fdrnew
     else
          fdrWPM_protective = ones(size(wpm_protective(1,:)));
     end

     %% compute FDR for WPMs (risk only)
     ind = find(wpm_pv_risk(1,:)<=0.05);

     if isempty(ind)~=1

          for i=1:length(ind)
               % combine local ranksum test and its p-value
               n1(i) = nnz(wpm_pv_risk(1,:)<=wpm_pv_risk(1,ind(i)) & wpm_risk(1,:)>=wpm_risk(1,ind(i)));
               n2(i) = nnz(wpm_pv_risk(2:end,:)<=wpm_pv_risk(1,ind(i)) & wpm_risk(2:end,:)>= wpm_risk(1,ind(i)))/N;
          end

          fdr = n2./n1;

          clear n1 n2

          fdrWPM_risk = ones(size(wpm_risk(1,:)));

          % sort fdr based on wpm_pv_risk and wpm_degree_local
          testM = [fdr' wpm_pv_risk(1,ind)' round(wpm_risk(1,ind)')];

          [testM idx] = sortrows(testM,-1);

          for i=1:size(testM,1)
               indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
               fdrnew(i) = min(testM(indtmp,1));
          end

          % assign FDR to WPMs
          fdrWPM_risk(ind(idx)) = fdrnew;
          clear fdr fdrnew
     else
          fdrWPM_risk = ones(size(wpm_risk(1,:)));
     end
     
     %% compute FDR for PATH (protective and risk)
     ind = find(path_pv(1,:)<=0.05);
     if isempty(ind)~=1
 
          for i=1:length(ind)
               n1(i) = nnz(path_pv(1,:)<=path_pv(1,ind(i)) & path(1,:)>=path(1,ind(i)));
               n2(i) = nnz(path_pv(2:end,:)<=path_pv(1,ind(i)) & path(2:end,:)>=path(1,ind(i)))/N; 
          end

          fdr = n2./n1;

          clear n1 n2
 
          fdrPATH = ones(size(path(1,:)));

          % sort fdr based on PATH scores
          testM = [fdr' path_pv(1,ind)' path(1,ind)'];
          [testM idx] = sortrows(testM,-1);

          for i=1:size(testM,1)
               indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
               fdrnew(i) = min(testM(indtmp,1));
          end

          % assign FDR to PATHs
          fdrPATH(ind(idx)) = fdrnew;

          clear fdr fdrnew
     else
          fdrPATH = ones(size(path(1,:)));
     end

     %% compute FDR for PATH (protective)
     ind = find(path_pv_protective(1,:)<=0.05);
     if isempty(ind)~=1

          for i=1:length(ind)
               n1(i) = nnz(path_pv_protective(1,:)<=path_pv_protective(1,ind(i)) & path_protective(1,:)>=path_protective(1,ind(i)));
               n2(i) = nnz(path_pv_protective(2:end,:)<=path_pv_protective(1,ind(i)) & path_protective(2:end,:)>=path_protective(1,ind(i)))/N;
          end

          fdr = n2./n1;

          clear n1 n2

          fdrPATH_protective = ones(size(path_protective(1,:)));

          % sort fdr based on PATH scores
          testM = [fdr' path_pv_protective(1,ind)' path_protective(1,ind)'];
          [testM idx] = sortrows(testM,-1);

          for i=1:size(testM,1)
               indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
               fdrnew(i) = min(testM(indtmp,1));
          end
          
          % assign FDR to PATHs
          fdrPATH_protective(ind(idx)) = fdrnew;
     
          clear fdr fdrnew
     else
          fdrPATH_protective = ones(size(path_protective(1,:)));
     end

      %% compute FDR for PATH (risk)
     ind = find(path_pv_risk(1,:)<=0.05);
     if isempty(ind)~=1

          for i=1:length(ind)
               n1(i) = nnz(path_pv_risk(1,:)<=path_pv_risk(1,ind(i)) & path_risk(1,:)>=path_risk(1,ind(i)));
               n2(i) = nnz(path_pv_risk(2:end,:)<=path_pv_risk(1,ind(i)) & path_risk(2:end,:)>=path_risk(1,ind(i)))/N;
          end

          fdr = n2./n1;

          clear n1 n2

          fdrPATH_risk = ones(size(path_risk(1,:)));

          % sort fdr based on PATH scores
          testM = [fdr' path_pv_risk(1,ind)' path_risk(1,ind)'];
          [testM idx] = sortrows(testM,-1);

          for i=1:size(testM,1)
               indtmp = find(testM(:,2)>=testM(i,2) & testM(:,3)<=testM(i,3));
               fdrnew(i) = min(testM(indtmp,1));
          end


          % assign FDR to PATHs
          fdrPATH_risk(ind(idx)) = fdrnew;
          clear fdr fdrnew

     else
          fdrPATH_risk = ones(size(path_risk(1,:)));
     end

     wpm_ranksum = wpm(1,:);
     path_ranksum = path(1,:);
     wpm_pv = wpm_pv(1,:);
     path_pv = path_pv(1,:);
end

bpm_ranksum = bpm(1,:);

bpm_pv = bpm_pv(1,:);

if exist('genesetname','var')==0
     outputfile = sprintf('results_%s_R0.mat',ssmFile);
else
     outputfile = sprintf('results_%s_%s_R0.mat',genesetname,ssmFile);
end

if exist('fdrBPM','var')|exist('fdrWPM','var')|exist('fdrPATH','var')
     save(outputfile,'fdr*','*_ranksum','*_pv')
end
