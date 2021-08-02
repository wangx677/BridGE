function computessi(model,marginal,alpha1,alpha2,plinkCluster2,nWorker,R)

% FUNCTION computessi(model,alpha1,alpha2,plinkCluster2,nWorker,R)
%
% COMPUTESSI computes SNP-SNP interactions based on given disease model and parameters.
% 
% INPUTS:
%   model - disease model assumptions
%   marginal - control individual SNP's marginal effect (joint mutation has to be more significant than the single SNP)
%              1 means control, 0 means no control
%   alpha1 - a significance constrain used in hygeSSI that controls joint mutation (11) significance 
%   alpha2 - a significance constrain used in hygeSSI that controls individual mutation (10,01) 
%          or wild type (00) signficance
%   plinkCluster2 - Output from plink size-2 clustering after removing individuals that
%           are not paired with others.
%   nWorker - number of workers for parallel computing 
%   R - permutation run ID. R = 0 is for real-phenotype data, R = 1,2,3,... are for randomized-phenotype data
% 
% OUTPUTS:
%  A mat-file include a matrix vector (size of 2) ssM. It's named based on given disease model. 
% 

if strcmp(model,'RR')
     if marginal==1
          ssmFile = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_RR',num2str(alpha1),num2str(alpha2));
     elseif marginal==0
          ssmFile = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_RR',num2str(alpha1),num2str(alpha2));
     end
elseif strcmp(model,'DD')
     if marginal==1
          ssmFile = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_DD',num2str(alpha1),num2str(alpha2));
     elseif marginal==0
          ssmFile = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_DD',num2str(alpha1),num2str(alpha2));
     end
elseif strcmp(model,'RD')
     if marginal==1
          ssmFile = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_RD',num2str(alpha1),num2str(alpha2));
     elseif marginal==0
          ssmFile = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_RD',num2str(alpha1),num2str(alpha2));
     end
elseif strcmp(model,'combined')
     if marginal==1
          ssmFile = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_combined',num2str(alpha1),num2str(alpha2));
     elseif marginal==0
          ssmFile = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_combined',num2str(alpha1),num2str(alpha2));
     end
elseif strcmp(model,'AA')
     ssmFile = 'ssM_lr_cassi_pv0.25';
end

% check if result file exist
ssmFileNew = sprintf('%s_R%s.mat',ssmFile,num2str(R));

if exist(ssmFileNew,'file')==2
     return;
end

% if ssmFileNew file doesn't exist
if nWorker==1
     kk = 1;
else
     dataFile='SNPdataAR.mat';
     load(dataFile)
     if length(SNPdata.rsid)>5000
          kk = 10; %define how to partition the data to invoke parallelized hygeSSI
     else
          kk = 5;
     end
end

% create parallel pool
pr = gcp('nocreate');
if isempty(pr)==1
     if nWorker>1
          parpool(nWorker)
          pr = 0;
     end
else
     pr = 1;
end

if R==0
     if strcmp(model,'AA')~=1
          parallelhygssi(model,marginal,alpha1,alpha2,nWorker,kk,R);
     elseif strcmp(model,'AA')
          % compute SNP-SNP interaction based on logistic regression
          system(sprintf('%s/scripts/runcassi.sh gwas_data_final LR 0.25 %s',getenv('BRIDGEPATH'),num2str(R)));
     end
else
     % ssM matrix from random phenotype labels
     dataFile='SNPdataAR.mat';
     load(dataFile)
     if isfield('SNPdata','randpheno')
          % use pre-defined random phenotype labels
          if isequal(sort(SNPdata.randpheno(:,R)),1:length(SNPdata.pheno))
               phenonew = SNPdata.randpheno(:,R);
          elseif isequal(sort(SNPdata.randpheno(R,:)),1:length(SNPdata.pheno))
               phenonew = SNPdata.randpheno(R,:);
          else
               error('Please check SNPdata.randpheno.')
          end
     else
          % if no pre-defined random phenotype label exists
          if isfield(SNPdata,'rand_pheno')
               % sample randomization is done indpendently
               phenonew = SNPdata.rand_pheno(:,R);
          elseif (nnz(SNPdata.pheno==0)~=nnz(SNPdata.pheno==1) & exist(plinkCluster2,'file')~=2)
               % data is not matched case-control
                rand('seed',R+1);
                phenonew = SNPdata.pheno(randperm(length(SNPdata.pheno)));
           elseif exist(plinkCluster2,'file')~=2
               rand('seed',R+1);
               phenonew = SNPdata.pheno(randperm(length(SNPdata.pheno)));
           else
                % randomized labels based on matched case-control shuffling
                phenonew = withinclassrand(dataFile,plinkCluster2,R+1);
           end
     end
            
     % compute ssM matrix with random phenotype labels 
     if strcmp(model,'AA')~=1
          parallelhygssi(model,marginal,alpha1,alpha2,nWorker,kk,R,phenonew);	
     elseif strcmp(model,'AA')
          % generate random plink file based on phenonew
          phenonew = phenonew+1; % in plink 1=control and 2=case
          FID = SNPdata.fid;
          PID = SNPdata.pid;
          if size(phenonew,1)==1
               phenonew = phenonew';
          end
          outputtmp = table(FID,PID,phenonew);
          writetable(outputtmp,'phenonewtmp.txt','Delimiter',' ','WriteVariableNames',0);
          system('plink --bfile gwas_data_final --noweb --pheno phenonewtmp.txt --make-bed --out plinktmp');

          % compute SNP-SNP interaction based on logistic regression               
          system(sprintf('%s/scripts/runcassi.sh %s LR 0.25 %s',getenv('BRIDGEPATH'),'plinktmp',num2str(R)));

          system('rm plinktmp.* phenonewtmp.txt')
     end
end

if pr==0
     % delete parallel pool if it is created by this function
     delete(gcp('nocreate'))
end

