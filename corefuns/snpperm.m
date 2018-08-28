function snpperm(ssmFile,bpmindFile,snpPerms,permBlockSize,minPath,nWorker,runNumber)

% FUNCTION snpperm(ssmFile,bpmindFile,snpPerms,permBlockSize,minPath,netDensity,ssMdensity,nWorker,runNumber)
%
% SNPPERM run SNP permutations for a sample population with fixed phenotype lables
%
% INPUTS:
%   ssmFile: SNP-SNP interaction matrix file
%   bpmindFile: BPM index file
%   snpPerms: total number of permutations
%   permBlockSize: number of permutations in each run
%   nWorker: number of matlab workers in parallel computing
%   runNumber:before total number of permutation is complete, an intermediate result will be saved 
%     in case there is a system failure. If that happens, runNumber tells from where to start again. 
%
% OUTPUTS:
%   A mat-file named as results_<ssmFile>_snpP_run<snpPerms/permBlockSize>_total<snpPerms>.mat 
% 

if (nargin > 9)
    sprintf('Incorrect number of input arguments!');
end


% default runNumber = 1
if exist('runNumber','var')~=1
    runNumber = 1;
end

% default nWorker = Maximum Number of Local Workers - 2
if exist('nWorker','var')~=1
    myCluster = parcluster('local');
    nWorker = myCluster.NumWorkers-2;
end

% if parallel pool not exist, create it
pr = gcp('nocreate');
if isempty(pr)==1
     parpool(nWorker)
     pr = 0;
else
     pr = 1;
end

% permBlockSize can't be greater than snpPerms
if permBlockSize>snpPerms
     permBlockSize = snpPerms;
end

% define the result file based on the SNP-SNP interaction matrix, netDensity, and total number of permutations 
resultFile = sprintf('results_%s_snpP_density%s_run%s_total%s.mat',ssmFile,num2str(netDensity),num2str(snpPerms/permBlockSize),num2str(snpPerms));

% if resultFile exists, load the file and return
if (exist(resultFile, 'file') == 2)
    load(resultFile)
    return;
end

% if resultFile doesn't exist, load SNP-SNP interaction matrix file
load(ssmFile);

% if SNP-SNP interaction matrix (ssM) is stored in vector format, convert it to matrix format
if (size(ssM{1},1)==1)
    for tt=1:2
        ssM{tt} = squareform(ssM{tt});
    end
end

% get the size of ssM
[p q] = size(ssM{1});

% convert sparse matrix to full matrix
for tt=1:2
    ssM{tt} = full(ssM{tt});
end

% compute chi-square stats for all BPMs/WPMs and ranksum stats for PATHs
genstatsFile=sprintf('genstats_%s.mat',ssmFile);

if (exist(genstatsFile, 'file') == 2)
    load(genstatsFile)
else
    [BPMranksum WPMranksum PATHranksum] = genstats(ssmFile,ssM,bpmindFile,minPath,nWorkers); 
end

% permutation initialization
if runNumber==1
     % if resultFile exist, load file
     tn=(runNumber)*permBlockSize;
     resultFile = sprintf('results_%s_snpP_density%s_run%s_total%s.mat',ssmFile,num2str(netDensity),num2str(runNumber),num2str(tn));
     if (exist(resultFile,'file') ==2)
          load(resultFile)
          load(bpmindFile)

          % remove fields that are not useful
          rmfield(BPM,{'path1idx','path2idx','ind1size','ind2size'})
          rmfield(WPM,'indsize');
     else
     
          % initial calclation
          load(bpmindFile); 

          % remove fields that are not useful
          rmfield(BPM,{'path1idx','path2idx','ind1size','ind2size'})
          rmfield(WPM,'indsize');

          if (size(ssM{1},1)==1)
               for tt=1:2
                    ssM{tt} = squareform(ssM{tt});
               end
          end
    
          % launch 1000 permutation for all BPM/WPM/PATH to estimate the empirical p-values before filtering.
          % although this step takes extra time, it has at least two purposes: 1) for validation; 2) for QQ plot 
  
          for tt = 1:2
               BPMpvalue{tt} = ones(1,length(BPM.size));
               WPMpvalue{tt} = ones(1,length(WPM.size));
               PATHpvalue{tt} = WPMpvalue{tt};
	       
               tmppermBlockSize = 1000;
               tmpIND2keep = 1:length(BPM.size);
               sumMM = sum(ssM{tt}); % for pathway ranksum test 	
               [BPMCounts{tt} WPMCounts{tt} PATHCounts{tt}] = ...
                        runsnppermutation(ssM{tt},tmppermBlockSize,BPM.size,WPM.size,p, ...
                            BPM.ind1,BPM.ind2,WPM.ind,WPM.ind,WPM.ind, ...
                        PATHranksum{tt},sumMM,tmpIND2keep,(runNumber-1)*permBlockSize);
              BPMpvalue{tt} = (BPMCounts{tt}+1)/tmppermBlockSize;
              WPMpvalue{tt} = (WPMCounts{tt}+1)/tmppermBlockSize;
              PATHpvalue{tt} = (PATHCounts{tt}+1)/tmppermBlockSize;
          end

          % filter BPMs/WPMs/PATHs (require marginal significance)
          cutoff = -log10(0.05);
          for tt=1:2
               IND2keepBPM{tt} = find(BPMchi2Local{tt}>=cutoff & BPMchi2{tt}>=cutoff & BPMpvalue{tt}'<=0.1);
               IND2keepWPM{tt} = find(WPMchi2Local{tt}>=cutoff & WPMchi2{tt}>=cutoff & WPMpvalue{tt}'<=0.1);
               IND2keepPATH{tt} = find(PATHranksum{tt}>=cutoff & PATHpvalue{tt}<=0.1); 
          end

          [p q]=size(ssM{1});

          % update information so that only do deep permutations for BPMs(in IND2keepBPM), WPM(inIND2keepWPM) PATHs(in IND2keepPATH)
          for tt=1:2
               % BPM
               BPMind1New{tt} = BPM.ind1(IND2keepBPM{tt});
               BPMind2New{tt} = BPM.ind2(IND2keepBPM{tt});
               BPMSizeNew{tt} = BPM.size(IND2keepBPM{tt});
               BPMdensityNew{tt} = BPMdensity{tt}(IND2keepBPM{tt});
               BPMchi2LocalNew{tt} = BPMchi2Local{tt}(IND2keepBPM{tt});
               SelfcomparisonCountsBPM{tt} = zeros(size(IND2keepBPM{tt}));
     
               % WPM
               WPMind = WPM.ind;
               WPMindNew{tt} = WPM.ind(IND2keepWPM{tt});
               WPMSizeNew{tt} = WPM.size(IND2keepWPM{tt});
               WPMdensityNew{tt} = WPMdensity{tt}(IND2keepWPM{tt});
               WPMchi2LocalNew{tt} = WPMchi2Local{tt}(IND2keepWPM{tt});
               SelfcomparisonCountsWPM{tt} = zeros(size(IND2keepWPM{tt}));

               % PATH
               PATHindNew{tt} = WPM.ind(IND2keepPATH{tt});
               PATHranksumNew{tt} = PATHranksum{tt}(IND2keepPATH{tt});
               SelfcomparisonCountsPATH{tt} = zeros(size(IND2keepPATH{tt}));
          end
          
          clear BPM WPM PATH BPMdensity BPMchi2Local WPMdensity WPMchi2Local PATHranksum tmp* 
     end
else
    % if file exists, load file
     tn=(runNumber-1)*permBlockSize;
     load(bpmindFile)
     load(sprintf('results_%s_snpP_density%s_run%s_total%s.mat',ssmFile,num2str(netDensity),num2str(runNumber-1),num2str(tn)))
     [p q]=size(ssM{1});
end

% if no BPM/WPM/PATH left after filtering, exit the program
if (length(IND2keepBPM{1})==0 & length(IND2keepBPM{2})==0 & length(IND2keepWPM{1})==0 & length(IND2keepWPM{2})==0 & length(IND2keepPATH{1})==0 & length(IND2keepPATH{2})==0);
     return;
end

% launch permutation
if permBlockSize>0
    while  runNumber<=snpPerms/permBlockSize
        totaln = runNumber*permBlockSize;
        resultFile = sprintf('results_%s_snpP_density%s_run%s_total%s.mat',ssmFile,num2str(netDensity),num2str(runNumber),num2str(totaln));
        
        if (exist(resultFile,'file') ==2)
            % if resultFile exists, load file
            load(resultFile)
            load(bpmindFile)
        else
            % if resultFile doesn't exists, run permutation
            for tt=1:2
                sumMM = sum(ssM{tt}); % for pathway ranksum test
                [BPMCounts{tt} WPMCounts{tt} PATHCounts{tt}] = ...
                        runsnppermutation(ssM{tt},permBlockSize,BPMSizeNew{tt},WPMSizeNew{tt},p, ...
                            BPMind1New{tt},BPMind2New{tt},WPMind,WPMindNew{tt},PATHindNew{tt},BPMdensityNew{tt}, ...
                      WPMdensityNew{tt},BPMchi2LocalNew{tt}, WPMchi2LocalNew{tt},...
                        PATHranksumNew{tt},sumMM,IND2keepBPM{tt},(runNumber-1)*permBlockSize);
                
                % update counts (random runs have better scores than real runs)
                SelfcomparisonCountsBPM{tt} =  SelfcomparisonCountsBPM{tt}+BPMCounts{tt}';
                SelfcomparisonCountsWPM{tt} = SelfcomparisonCountsWPM{tt}+WPMCounts{tt}';
                SelfcomparisonCountsPATH{tt} = SelfcomparisonCountsPATH{tt}+PATHCounts{tt};

                % update permutation p-values
                BPMpvalue{tt}(IND2keepBPM{tt}) = (SelfcomparisonCountsBPM{tt}+1)/totaln;
                WPMpvalue{tt}(IND2keepWPM{tt}) = (SelfcomparisonCountsWPM{tt}+1)/totaln;
                PATHpvalue{tt}(IND2keepPATH{tt}) = (SelfcomparisonCountsPATH{tt}+1)/totaln;
            end
         
            % save information to resultFile
            save(sprintf('results_%s_snpP_density%s_run%s_total%s.mat',ssmFile,num2str(netDensity),num2str(runNumber),num2str(totaln)), ...
               '*pvalue','ssM','ssMdensity','ssMnetcut*','Selfcomparison*','IND2keep*','*New','WPMind','-v7.3')

            % delete intermediate resultFiles
            delete(sprintf('results_%s_snpP_density%s_run%s_total%s.mat',ssmFile,num2str(netDensity),num2str(runNumber-2),num2str((runNumber-2)*permBlockSize)))
        end
        runNumber = runNumber+1;
    end
end

if pr==0
     % delete parallel pool if it is created by this function
     delete(gcp('nocreate'))
end

end
