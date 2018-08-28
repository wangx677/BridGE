function ssM = hygessiwithin(data,pheno,alpha1,alpha2)

% FUNCTION ssM = hygessiwithin(data,pheno,alpha1,alpha2)
% 
% HYGESSIWITHIN compute pairwise using same set of data (e.g. recessive-recessive 
% or dominant-dominant SNP-SNP interactions)
%
% INPUTS:
%   datar - SNP data in recessive format 
%   datad - SNP data in dominant format
%   pheno - phenotype labels
%   alpha1 - a signifcance constrain used in hygeSSI that controls joint mutation (11) significance
%   alpha2 - a signifcance constrain used in hygeSSI that controls individual mutation (10,01)
%      or wild type (00) signficance
% 
% OUTPUTS:
%   A mat-file include a matrix vector (size of 2) ssM. ssM{1} is for protective interactions and 
%      ssM{2} is for risk interactions.
%

if (nargin > 4)
    sprintf('Incorrect number of input arguments!');
end

if exist('hyge_dictionary_table.mat','file')==2
     load hyge_dictionary_table.mat
else
     HygeTable = hygedictionary(length(pheno),nnz(pheno==1),nnz(pheno==0));
end

x = find(pheno==0); % controls
y = find(pheno==1); % cases

% setup data
data1 = data;
data0 = 1- data;

case1 = data(y,:);
case0 = 1-case1;

control1 = data(x,:);
control0 = 1-control1;

[p q] = size(data);

clear data

tmp = ones(q,q);
for tt=1:2
     ssM{tt} = sparse(q,q);
     pair2keep{tt} = triu(tmp,1)+tril(tmp,-1); % diagonal all zeros
end

clear tmp

%% compute individual SNP effect
% recessive disease model
hyge1Protective = zeros(1,q);
hyge1Risk = zeros(1,q);

a = sum(data1(x,:));
b = sum(data1(y,:));

indtmp = find(a+b>0);

test = [ones(length(indtmp),1)*length(x), a(indtmp)'+b(indtmp)',a(indtmp)'];
[tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
hyge1Protective(indtmp) = HygeTable(ind,4);

test = [ones(length(indtmp),1)*length(y), a(indtmp)'+b(indtmp)',b(indtmp)'];
[tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
hyge1Risk(indtmp) = HygeTable(ind,4);

clear a b indtmp test tmp ind  

%% compute snp pair (11) effect

C = data1'*data1;
C = squareform(tril(C,-1) + triu(C,1));
Ctmp = C;
Ctmp(squareform(pair2keep{1})+squareform(pair2keep{2})==0) = [];
uC = unique(Ctmp);
uC(uC==0) = [];

[uCutx uCuty] = ucuthygeinv(x,y,alpha1,p,uC);

[testMx testMy] = replacehyge(x,y,C,uC,uCutx,uCuty,p);

Scase = case1'*case1;
Scontrol = control1'*control1;
Scase = tril(Scase,-1)+triu(Scase,1);
Scontrol = tril(Scontrol,-1)+triu(Scontrol,1);

pair2keep{1}(Scontrol-(squareform(testMx)+1)<0) = 0;
pair2keep{2}(Scase-(squareform(testMy)+1)<0) = 0;

clear testMx testMy C Ctmp uC uCutx uCuty Scase Scontrol

% check significance of "11","00","01","10"combination 

for tt=1:2
     tmp = nchoosek(q,2);
     ssM11{tt} = zeros(1,tmp);
     ssM00{tt} = zeros(1,tmp);
     ssM10{tt} = zeros(1,tmp);
     ssM01{tt} = zeros(1,tmp);
     clear tmp

     IND{tt} = find(squareform(pair2keep{tt})==1);

     a = data1(x,:)'*data1(x,:);
     b = data1(y,:)'*data1(y,:);
     a = triu(a,+1)+tril(a,-1);
     b = triu(b,+1)+tril(b,-1);
     a = squareform(a);
     b = squareform(b);
     

     if (tt==1)
          ssM11{tt} = computessm(ssM11{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM11{tt} = computessm(ssM11{tt},IND{tt},y,b,a,HygeTable);
     end     
     ssM11{tt} = sparse(ssM11{tt});
      
     a = data0(x,:)'*data0(x,:);
     b = data0(y,:)'*data0(y,:);
     a = triu(a,+1)+tril(a,-1);
     b = triu(b,+1)+tril(b,-1);
     a = squareform(a);
     b = squareform(b);
     if (tt==1)
          ssM00{tt} = computessm(ssM00{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM00{tt} = computessm(ssM00{tt},IND{tt},y,b,a,HygeTable);
     end
     ssM00{tt} = sparse(ssM00{tt});

     a = data1(x,:)'*data0(x,:);
     b = data1(y,:)'*data0(y,:);
     a = squareform(triu(a,1)'); 
     b = squareform(triu(b,1)');
     if (tt==1)
          ssM10{tt} = computessm(ssM10{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM10{tt} = computessm(ssM10{tt},IND{tt},y,b,a,HygeTable);
     end
     ssM10{tt} = sparse(ssM10{tt});
     
     a = data1(x,:)'*data0(x,:);
     b = data1(y,:)'*data0(y,:);
     a = squareform(tril(a,-1));
     b = squareform(tril(b,-1));
     if (tt==1)
          ssM01{tt} = computessm(ssM01{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM01{tt} = computessm(ssM01{tt},IND{tt},y,b,a,HygeTable);
     end
     ssM01{tt} = sparse(ssM01{tt});
end

clear a b IND HygeTable

%% apply hygeSSI

hyge1{1} = max(repmat(hyge1Protective',1,q),repmat(hyge1Protective,q,1));
hyge1{2} = max(repmat(hyge1Risk',1,q),repmat(hyge1Risk,q,1));

hyge1{1} = tril(hyge1{1},-1) + triu(hyge1{1},1);
hyge1{2} = tril(hyge1{2},-1) + triu(hyge1{2},1);

hyge1{1} = squareform(hyge1{1});
hyge1{2} = squareform(hyge1{2});

for tt=1:2
     ssM{tt} = (ssM11{tt} - max(ssM00{tt},max(ssM10{tt},ssM01{tt}))).*squareform(pair2keep{tt});
end

for tt=1:2
     % set negative scores to be 0
     ssM{tt}(ssM{tt}<0) = 0;
     % apply constrains
     ssM{tt}(ssM00{tt}>=-log10(alpha2)) = 0;
     ssM{tt}(ssM01{tt}>=-log10(alpha2)) = 0;
     ssM{tt}(ssM10{tt}>=-log10(alpha2)) = 0;
     % convert to matrix format
     ssM{tt} = squareform(ssM{tt});
end

end

function [testMx testMy] = replacehyge(x,y,C,uC,uCutx,uCuty,p)
% This function generates two matrices that replaces all uC values in C with 
% corresponding uCutx, or uCuty
testMx = replace(C,uC',uCutx');
testMx(testMx==0) = p;

if length(x) == length(y)
      testMy = testMx;
else
      testMy = replace(C,uC',uCuty');
      testMy(testMy==0) = p;
end
end


function [uCutx uCuty] = ucuthygeinv(x,y,alpha,p,uC)
% This function use inverse hypergeometric to calculate given total number 
% of samples p, total number of events uC, and an signficance level alpha, 
% what's the required number events in x or y. 
uCutx = hygeinv(1-alpha,p,uC,length(x));

if length(x) == length(y)
      uCuty = uCutx;
else
      uCuty = hygeinv(1-alpha,p,uC,length(y));
end
end

function ssM = computessm(ssM,IND,x,a,b,HygeTable)
% compute snp-snp interactions scores based on HygeTable
test = [ones(length(IND),1)*length(x) a(IND)'+b(IND)' a(IND)'];
[tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
if (nnz(ind==0)>0)
      ind1 = find(ind~=0 & test(:,1)~=0);
      ind = ind(ind1);
      indnew = IND(ind1);
      ssM(indnew) = HygeTable(ind,4);
else
      ssM(IND) = HygeTable(ind,4);
end
end

