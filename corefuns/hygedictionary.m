function HygeTable = hygedictionary(N,n1,n2,nWorker)

% FUNCTION HygeTable = hygedictionary(data,pheno)
% 
% HYGEDICTIONARY computes all possible hypergeometric statistics given 
% the sample phenotype distribution. The result table can be used as a 
% dictionary in later steps
%
% INPUTS:
% N - total number of samples
% n1 - total number of cases 
% n2 - total number of controls
% nWorker - number of workers for parallel computing (optional)
% 
% OUTPUTS:
% HygeTable - a matrix with 4 columns. The total number of samples is fixed, 
%    so only the following information are included in the HygeTable: 
%    1) number of samples selected
%    2) total number of samples with mutations
%    3) number of selected samples have 
%

if (nargin < 3 | nargin >4)
    sprintf('Incorrect number of input arguments');
end


uC = 1:N; % all possible number of samples with mutations


if exist('nWorker','var')
     % create parallel pool
     pr = gcp('nocreate');

     if isempty(pr) == 1
          parpool(nWorker)
          pr = 0;
     else
          pr = 1;
     end

     if (n1==n2)
          parfor i=1:length(uC)
               k = 0:uC(i);
               tmp = hygetest(N,n1,k,uC(i));
               tmp = [repmat(n1,length(k),1) repmat(uC(i),length(k),1) k' tmp'];
               tmpHygeTable{i} = tmp;
          end
     else
          parfor i=1:length(uC)
                k = 0:uC(i);
                tmp1 = hygetest(N,n1,k,uC(i));
                tmp1 = [repmat(n1,length(k),1) repmat(uC(i),length(k),1) k' tmp1'];
                tmp2 = hygetest(N,n2,k,uC(i));
                tmp2 = [repmat(n2,length(k),1) repmat(uC(i),length(k),1) k' tmp2'];
                tmpHygeTable{i} = [tmp1;tmp2];
         end
     end

     HygeTable = cell2mat(tmpHygeTable');

     if pr==0
          % delete parallel pool if it is created by this function
          delete(gcp('nocreate'))
     end

else
    % non-parallel 
    HygeTable = [];

    if (n1==n2)
       for i=1:length(uC)
          k = 0:uC(i);
          tmp = hygetest(N,n1,k,uC(i));
          tmp = [repmat(n1,length(k),1) repmat(uC(i),length(k),1) k' tmp'];
          HygeTable = [HygeTable;tmp];
       end
     else
       for i=1:length(uC)
                k = 0:uC(i);
                tmp1 = hygetest(N,n1,k,uC(i));
                tmp1 = [repmat(n1,length(k),1) repmat(uC(i),length(k),1) k' tmp1'];
                tmp2 = hygetest(N,n2,k,uC(i));
                tmp2 = [repmat(n2,length(k),1) repmat(uC(i),length(k),1) k' tmp2'];
                HygeTable = [HygeTable;tmp1;tmp2];
        end
     end
end 

save hyge_dictionary_table HygeTable -v7.3
