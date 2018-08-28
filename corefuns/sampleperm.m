function sampleperm(ssmFile,bpmindFile,netDensity,snpPerms,permBlockSize,minPath,nWorker,R)

% FUNCTION sampleperm(ssmFile,bpmindFile,netDensity,snpPerms,permBlockSize,minPath,nWorker,R)
%
% SAMPLEPERMUTATION run sample permutations with nested SNP permutations in each run
% 
% INPUTS:
%   ssmFile - SNP-SNP interaction file name without suffix '_R<R>.mat'
%   bpmindFile - a mat-file that includes SNP indexes of pathways for all possbile BPMs, WPMs
%   netDensity - overall network density (protective and risk interactions combined)
%   snpPerms - number of SNP permutations in each run
%   permBlockSize - permutation block size
%   minPath - minimum pathway size
%   nWorker - number of workers for parallel computing 
%   R - permutation run ID. R = 0 is for real-phenotype data, R = 1,2,3,... are for randomized-phenotype data
% 
% OUTPUTS:
%  An mat-file named results_<ssmFile>_R<R>.mat_snpP_density<netDensity>_run<snpPerms/permBlockSize>_total<snpPerms>.mat'

if (nargin > 13)
    sprintf('Incorrect number of input arguments!');
end

% permBlockSize can't be greater than snpPerms
if (permBlockSize>snpPerms)
     permBlockSize = snpPerms;
end

% check if result file exist
ssmFileNew = sprintf('%s_R%s.mat',ssmFile,num2str(R));
samplePermFile = sprintf('results_%s_snpP_density%s_run%s_total%s.mat',ssmFileNew,num2str(netDensity),num2str(snpPerms/permBlockSize),num2str(snpPerms));

if exist(samplePermFile,'file')==2
     return;
end

% if sample permutation file doesn't exist
% create parallel pool
pr = gcp('nocreate');
if isempty(pr)==1
     parpool(nWorker)
     pr = 0;
else
     pr = 1;
end

% netDensity is the density for protective/risk combined network, 
% we need to find out the individual density (ssMdensity) for corresponding
% protective and risk networks. 

% if genstats file exist, get ssMdenisty information by loading this file
tmpfile = sprintf('genstats_%s_density%s',ssmFile,num2str(netDensity));
if exist(tmpfile,'file')==2
     load(tmpfile)
else 
     % else load ssM matrix from real data and compute ssMdensity 
     ssmFile0 = sprintf('%s_R%s.mat',ssmFile,num2str(0));

     if exist(ssmFile0,'file')
          load(ssmFile0);
     elseif R==0
          sprintf('%s not exist!',ssmFile0);
          return;
     end

     if (size(ssM{1},1)==1)
          for tt=1:2
               ssM{tt} = squareform(ssM{tt});
          end
     end

     [p q] = size(ssM{1});
     tmp = max(ssM{1},ssM{2});
     netcutALL = quantile(tmp(:),1-netDensity);

     for tt=1:2
          ssM{tt} = ssM{tt}>=netcutALL;
          ssMdensity{tt} = nnz(ssM{tt})/(p*q);
     end
end
     
% only keep variables that are needed
clearvars -except pr ssMdensity dataFile ssmFile ssmFileNew bpmindFile netDensity snpPerms permBlockSize minPath nWorker R

% run SNP permutation 
if R==0
     % netDensity is the density for protective/risk combined network, 
     % we need to find out the individual density (densityReal) for protective and risk network
     % This is can be derived after exploring the real label data  
     if exist('ssMdensity','var')~=1
          ssMdensity=[];
     end
     snpperm(ssmFileNew,bpmindFile,snpPerms,permBlockSize,minPath,netDensity,ssMdensity,nWorker,1);

else
     if (exist('ssMdensity','var')~=1) | (isempty(ssMdensity)==1)
          sprintf('Missing ssMdensity for random run R%s!',num2str(R))
     else
          % for data with real phenotype labels, ssMdensity can be missing or empty
          % but for data with random phenotype labels, ssMdensity can't.
          % Use a different name to differentiate them
          ssMdensityReal = ssMdensity;
     end

     snpperm(ssmFileNew,bpmindFile,snpPerms,permBlockSize,minPath,netDensity,ssMdensityReal,nWorker,1);
end

if pr==0
     % delete parallel pool if it is created by this function
     delete(gcp('nocreate'))
end

