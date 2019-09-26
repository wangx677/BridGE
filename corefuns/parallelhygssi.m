function ssM = parallelhygssi(model,marginal,alpha1,alpha2,nWorker,nMatrices,R,NewPheno)

% FUNCTION ssM = parallelhygssi(model,alpha1,alpha2,nWorker,k,R,NewPheno)
%
% PARALLELHYGESSI computes all pairwise hygeSSI scores based on given disease 
%   model assumptions (a SNP is recessive or dominant). 
%
% INPUTS:
% model - disease model assumptions. 
%    RR: both SNPs are recessive,
%    DD: both SNPs are dominant;
%    RD: one SNP is recessive and another is dominant
%    combined: take maximum score of 1,2,3 as interaction score
%    marginal - control individual SNP's marginal effect (joint mutation has to be more significant than the single SNP)
%              1 means control, 0 means no control
% alpha1 - a signifcance constrain used in hygeSSI that controls joint mutation (11) significance
% alpha2 - a signifcance constrain used in hygeSSI that controls individual mutation (10,01) 
%          or wild type (00) signficance
% nWorker - number of matlab workers for parallel computing
% nMatrices - divide SNP data to N matrices for parallel computing
% R - an index that is used to label ssM matrix file name. 0 for real pheno labels and others for
%     random pheno labels
% NewPheno - randmized pheno labels
%
% OUTPUTS:
% A mat-file include a matrix vector (size of 2) ssM. ssM{1} is for protective interactions and 
% ssM{2} is for risk interactions.The mat-file is named based on give disease model assumptions (model):
%   1, recessive-recessive, ssM_hygeSSI_alpha1<alpha1>_alpha2<alpha2>_RR_R<R>.mat
%   2, dominant-dominant, ssM_hygeSSI_alpha1<alpha1>_alpha2<alpha2>_DD_R<R>.mat
%   3, recessive-dominant, ssM_hygeSSI_alpha1<alpha1>_alpha2<alpha2>_RD_R<R>.mat
%   4, combined, ssM_hygeSSI_alpha1<alpha1>_alpha2<alpha2>_combined_R<R>.mat
%

if (nargin > 7)
    sprintf('Incorrect number of input arguments!');
end

if exist('R','var')==0
     R = 0;
end

pr = gcp('nocreate');
if isempty(pr)==1
     % if parallel pool is empty, create pool 
     if nWorker>1
          parpool(nWorker)
          pr = 0;
     end
else
     % if parallel pool exists, use current pool
     pr = 1; 
end

% define file name
switch model
case 'RR' % recessive-recessive interaction
     if marginal==0
          filename = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_RR_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     elseif marginal==1
          filename = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_RR_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     end
case 'DD' % dominant-dominant interaction
     if marginal==0
          filename = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_DD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     elseif marginal==1
          filename = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_DD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     end
case 'RD' % recessive-dominant interaction
     if marginal==0
          filename = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_RD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     elseif marginal==1
          filename = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_RD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     end
case 'combined' % combined interaction
     if marginal==0
          filename = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_combined_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     elseif marginal==1
          filename = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_combined_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
     end
end

if exist(filename,'file')==0
     load SNPdataAD.mat; % load a SNPdata file to get the number of SNPs

     nSNP = size(SNPdata.data,2);
     clear SNPdata

     % partition data for parallelization
     % if number of SNPs is not a multiple of nMatrices
     if nSNP/nMatrices ~= floor(nSNP/nMatrices) 
          idx = 1:floor(nSNP/nMatrices):nSNP;
          ind1 = idx(1:length(idx)-1);
          ind2 = idx(2:length(idx))-1;
          ind2(end) = nSNP;
     else
     % if number of SNPs is a multiple of nMatrices
          idx = 1:nSNP/nMatrices:nSNP;
          ind1 = idx(1:length(idx));
          ind2 = idx(2:length(idx))-1;
          ind2(length(ind2)+1) = nSNP;
     end

     if strcmp(model,'RD')
          rows1 = repmat(ind1',1,nMatrices);
          rows2 = repmat(ind2',1,nMatrices);
          cols1 = repmat(ind1,nMatrices,1);
          cols2 = repmat(ind2,nMatrices,1);
     else
          rows1 = triu(repmat(ind1',1,nMatrices));
          rows2 = triu(repmat(ind2',1,nMatrices));
          cols1 = triu(repmat(ind1,nMatrices,1));
          cols2 = triu(repmat(ind2,nMatrices,1));
     end

     I = find(rows1~=0);
     i1 = rows1(I);
     i2 = rows2(I);
     j1 = cols1(I);
     j2 = cols2(I);

     N = length(I);

     clear nSNP
     
     switch model

     case 'RR'  
          % load data
          load SNPdataAR.mat
          data = SNPdata.data;
          pheno = SNPdata.pheno;
          clear SNPdata
          
          % change the pheno if there is new pheno available         
          if exist('NewPheno','var')
               pheno = NewPheno;
          end
          
          % parallel computing
          ssM = parforhygessi(data,data,i1,i2,j1,j2,pheno,marginal,alpha1,alpha2,N);
          save(filename,'ssM','-v7.3');

     case 'DD'
          % load data
          load SNPdataAD.mat;
          data = SNPdata.data;
          pheno = SNPdata.pheno;
          clear SNPdata
          
          % change the pheno if there is new pheno available          
          if exist('NewPheno','var')
               pheno = NewPheno;
          end
          
          % parallel computing
          ssM = parforhygessi(data,data,i1,i2,j1,j2,pheno,marginal,alpha1,alpha2,N);
          save(filename,'ssM','-v7.3');

     case 'RD'
          % load data
          load SNPdataAR.mat; datar = SNPdata.data;
          load SNPdataAD.mat; datad = SNPdata.data;
          pheno = SNPdata.pheno;

          % change the pheno if there is new pheno available          
          if exist('NewPheno','var')
               pheno = NewPheno;
          end

          % parallel computing
          ssM = parforhygessi(datar,datad,i1,i2,j1,j2,pheno,marginal,alpha1,alpha2,N);
          save(filename,'ssM','-v7.3');

     case 'combined'
          if exist('NewPheno','var')
               pheno = NewPheno;
          else
               load SNPdataAR.mat
               pheno = SNPdata.pheno;
          end
          
          % if recessive-recessive, dominant-dominant, recessive-dominant interaction files don't exist,
          % call parallelhygssi to compute them
          if marginal==1
               file1 = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_RR_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
          elseif marginal==0
               file1 = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_RR_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
          end
          if exist(file1,'file')==0
               ssM1 = parallelhygssi('RR',marginal,alpha1,alpha2,nWorker,nMatrices,R,pheno);
          else
               load(file1); ssM1 = ssM; clear ssM
          end

          if marginal==1
               file2 = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_DD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
          elseif marginal==0
               file2 = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_DD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
          end

          if exist(file2,'file')==0
               ssM2 = parallelhygssi('DD',marginal,alpha1,alpha2,nWorker,nMatrices,R,pheno);
          else
               load(file2); ssM2 = ssM; clear ssM
          end

          if marginal==1
               file3 = sprintf('ssM_hygeSSI_alpha1%s_alpha2%s_RD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
          elseif marginal==0
               file3 = sprintf('ssM_mhygeSSI_alpha1%s_alpha2%s_RD_R%s.mat',num2str(alpha1),num2str(alpha2),num2str(R));
          end
          if exist(file3,'file')==0
               ssM3 = parallelhygssi('RD',marginal,alpha1,alpha2,nWorker,nMatrices,R,pheno);              
          else
		       load(file3); ssM3 = ssM; clear ssM
          end
          
          % combine three matrices
          [ssM maxidx] = ssicombine(ssM1,ssM2,ssM3)
          
          save(filename,'ssM','maxidx','-v7.3')
     end
else
     load(filename);
end 

if pr==0
     % delete parallel pool if it is created by this function
     delete(gcp('nocreate'))
end

end


function ssM = parforhygessi(datar,datad,i1,i2,j1,j2,pheno,marginal,alpha1,alpha2,N)

% FUNCTION parforhygessi(datar,datad,i1,i2,j1,j2,pheno,marginal,alpha1,alpha2,N)
%
% PARFORHYGESSI use parallelized loops to compute pairwise interactions in sub-matrices

parfor i=1:N
	ssMtmp{i} = hygessiacross(datar(:,i1(i):i2(i)),datad(:,j1(i):j2(i)),pheno,marginal,alpha1,alpha2);
end

for i=1:N
	for tt=1:2
        	ssM{tt}(i1(i):i2(i),j1(i):j2(i)) = ssMtmp{i}{tt};
        end
end

for tt=1:2
	ssM{tt} = max(ssM{tt},ssM{tt}');
	ssM{tt} = triu(ssM{tt},1)+tril(ssM{tt},-1);
        ssM{tt} = sparse(squareform(ssM{tt}));
end
end


function [ssM maxidx] = ssicombine(ssM1,ssM2,ssM3)

% FUNCTION [ssM maxidx] = ssicombine(ssM1,ssM2,ssM3,output)
% 
% SSICOMBINE combines three interaction matrices and keep track of which disease model
% gives the maximum interaction score for each pair

for tt=1:2
        maxidx{tt} = zeros(size(ssM1{tt}));
        ssMtmp{tt} = max(ssM1{tt},ssM2{tt});
        maxidx{tt}(ssM1{tt} >= ssM2{tt}) = 1;
        maxidx{tt}(ssM2{tt} > ssM1{tt}) = 2; 
        ssM{tt} = max(ssMtmp{tt},ssM3{tt});
        maxidx{tt}(ssM3{tt} > ssMtmp{tt}) = 3;
        maxidx{tt} = maxidx{tt}.*(ssM{tt}>0);
end

end
