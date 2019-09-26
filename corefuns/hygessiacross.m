function ssM = hygessiacross(datar,datad,pheno,marginal,alpha1,alpha2)

% FUNCTION ssM = hygessiacross(datar,datad,pheno,marginal,alpha1,alpha2)
% 
% HYGESSIACROSS compute pairwise interactions between two datasets (e.g. recessive and dominant SNPs)
%
% IMPORTANT NOTICE: the final ssM matrix should be symmetric, however when it's 
% implemented for parallel computing, the imput data is just a subset of the 
% original data. So we didn't adjust the ssM matrix to make it symmetric.
%
% INPUTS:
%   datar - SNP data in recessive format 
%   datad - SNP data in dominant format
%   pheno - phenotype labels
%   marginal - control individual SNP's marginal effect (joint mutation has to be more significant than the single SNP)
%              1 means control, 0 means no control
%   alpha1 - a signifcance constrain used in hygeSSI that controls joint mutation (11) significance
%   alpha2 - a signifcance constrain used in hygeSSI that controls individual mutation (10,01)
%     or wild type (00) signficance
% 
% OUTPUTS:
%   A mat-file include a matrix vector (size of 2) ssM. ssM{1} is for protective interactions and 
%     ssM{2} is for risk interactions.
%

if (nargin > 5)
    sprintf('Incorrect number of input arguments!');
end

% if datar and datad are identical, run hygessiwithin 
if isequal(datar,datad)==1
    ssM = hygessiwithin(datar,pheno,marginal,alpha1,alpha2);
    return;
end

% check if hypergeometric dictionary fie exists,
% if not, compute all possible hypergeometric stats for given sample size

if exist('hyge_dictionary_table.mat','file')==2
     load('hyge_dictionary_table.mat')
else
     HygeTable = hygedictionary(length(pheno),nnz(pheno==1),nnz(pheno==0));
end

x = find(pheno==0);
y = find(pheno==1);

% setup data
datar1 = datar;
datar0 = 1- datar;

caser1 = datar1(y,:);
controlr1 = datar1(x,:);

datad1 = datad;
datad0 = 1- datad;

cased1 = datad1(y,:);
controld1 = datad1(x,:);

p = size(datar,2);
q = size(datad,2);
n = size(datar,1);

clear datar datad

for tt=1:2
     ssM{tt} = sparse(p,q); % snp-snp interaction matrices
     pair2keep{tt} = ones(p,q); % snp-snp pairs to keep
end

%% compute individual SNP effect
% recessive disease model
hyge1ProtectiveRecessive = zeros(1,p); % mutations are more frequent in controls
hyge1RiskRecessive = zeros(1,p); % mutations are more frequent in cases

a = sum(datar1(x,:));
b = sum(datar1(y,:));

% find SNPs that with at least 1 nonzeros
indtmp = find(a+b>0);

% mutations with protecitve effect
% create a three column matrix: number of selected samples, total number of 
% infected samples, number of infected samples in selected samples 
if length(indtmp)>0
     test = [ones(length(indtmp),1)*length(x), a(indtmp)'+b(indtmp)',a(indtmp)'];
     % find the index of these combinations in the dictionay table
     [tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
     % assign the hypergeometric stats to given data
     hyge1ProtectiveRecessive(indtmp) = HygeTable(ind,4);
end

% mutations with risk effect, 
if length(indtmp)>0
     test = [ones(length(indtmp),1)*length(y), a(indtmp)'+b(indtmp)',b(indtmp)'];
     [tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
     hyge1RiskRecessive(indtmp) = HygeTable(ind,4);
end

% dominant disease model 
hyge1ProtectiveDominant = zeros(1,q);
hyge1RiskDominant = zeros(1,q);

a = sum(datad1(x,:));
b = sum(datad1(y,:));

indtmp = find(a+b>0);

if length(indtmp)>0
     test = [ones(length(indtmp),1)*length(x), a(indtmp)'+b(indtmp)',a(indtmp)'];
     [tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
     hyge1ProtectiveDominant(indtmp) = HygeTable(ind,4);
end

if length(indtmp)>0
     test = [ones(length(indtmp),1)*length(y), a(indtmp)'+b(indtmp)',b(indtmp)'];
     [tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
     hyge1RiskDominant(indtmp) = HygeTable(ind,4);
end

clear a b test tmp ind indtmp

%% compute snp pair (11) effect

% 11 frequencies in all samples
C = datar1'*datad1;
% unique 11 frequencies (no need to compute for duplicates)
uC = unique(C);
% no need to compute when 11 frequency = 0
uC(uC==0) = [];
% use reverse hypergeometric function to find out, to achieve the signficance cutoff alpha1, 
% how many selected samples need to be infected (11)
[uCutx, uCuty] = ucuthygeinv(x,y,alpha1,n,uC);

% construct test matrices that reflect number of selected samples 
% need to be infected (11) to reach the signficance cutoff alpha1
[testMx, testMy] = replacehyge(x,y,C,uC,uCutx,uCuty,n);

% number of cases have 11
Scase = caser1'*cased1;
% number of controls have 11
Scontrol=controlr1'*controld1;

% only keep SNP pair whoes the number of observed infected samples (11) in 
% controls/cases is greater than the required numbers, 
pair2keep{1}(Scontrol-(testMx+1)<0) = 0;
pair2keep{2}(Scase-(testMy+1)<0) = 0;

clear C uC uC uCutx uCuty testMx testMy Scase Scontrol

% check significance of "11","00","01","10" combination 
for tt=1:2
     ssM11{tt} = zeros(p*q,1);
     ssM00{tt} = zeros(p*q,1);
     ssM01{tt} = zeros(p*q,1);
     ssM10{tt} = zeros(p*q,1);

     pair2keep{tt} = pair2keep{tt}(:);

     IND{tt} = find(pair2keep{tt}==1);

     a = datar1(x,:)'*datad1(x,:);
     b = datar1(y,:)'*datad1(y,:);
     if (tt==1)
          ssM11{tt} = computessm(ssM11{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM11{tt} = computessm(ssM11{tt},IND{tt},y,b,a,HygeTable);
     end    
     ssM11{tt} = sparse(reshape(ssM11{tt},p,q));
            
     a = datar0(x,:)'*datad0(x,:);
     b = datar0(y,:)'*datad0(y,:);          
     if (tt==1)
          ssM00{tt} = computessm(ssM00{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM00{tt} = computessm(ssM00{tt},IND{tt},y,b,a,HygeTable);
     end
     ssM00{tt} = sparse(reshape(ssM00{tt},p,q));

     a = datar1(x,:)'*datad0(x,:);
     b = datar1(y,:)'*datad0(y,:);
     if (tt==1)
          ssM10{tt} = computessm(ssM10{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM10{tt} = computessm(ssM10{tt},IND{tt},y,b,a,HygeTable);
     end     
     ssM01{tt} = sparse(reshape(ssM01{tt},p,q));

     a = datar0(x,:)'*datad1(x,:);
     b = datar0(y,:)'*datad1(y,:);
     if (tt==1)
          ssM01{tt} = computessm(ssM01{tt},IND{tt},x,a,b,HygeTable);
     else
          ssM01{tt} = computessm(ssM01{tt},IND{tt},y,b,a,HygeTable);
     end
     ssM10{tt} = sparse(reshape(ssM10{tt},p,q));

     pair2keep{tt} = reshape(pair2keep{tt},p,q);
end

clear a b IND HygeTable

%% apply hygeSSI
% individual SNP effect
hyge1{1} = max(repmat(hyge1ProtectiveRecessive',1,q),repmat(hyge1ProtectiveDominant,p,1));
hyge1{2} = max(repmat(hyge1RiskRecessive',1,q),repmat(hyge1RiskDominant,p,1));

% ssM11 - max(ssM00, ssM01, ssM10, hyge1)
for tt=1:2
     ssM00{tt}(isinf(ssM00{tt})) = 16;
     ssM11{tt}(isinf(ssM11{tt})) = 16;
     ssM10{tt}(isinf(ssM10{tt})) = 16;
     ssM01{tt}(isinf(ssM01{tt})) = 16;
     hyge1{tt}(isinf(hyge1{tt})) = 16;
     
     if marginal==0
          ssM{tt} = (ssM11{tt}-max(ssM00{tt},max(ssM01{tt},ssM10{tt}))).*pair2keep{tt};
     elseif marginal==1
          ssM{tt} = (ssM11{tt}-max(ssM00{tt},max(ssM01{tt},max(ssM10{tt},hyge1{tt})))).*pair2keep{tt};
     end
end
 
for tt=1:2
     % set negative scores to be 0
     ssM{tt}(ssM{tt}<0) = 0; 
     % apply constrains
     ssM{tt}(ssM00{tt}>=-log10(alpha2)) = 0; 
     ssM{tt}(ssM01{tt}>=-log10(alpha2)) = 0;
     ssM{tt}(ssM10{tt}>=-log10(alpha2)) = 0;

     % IMPORTANT NOTICE
     % the final ssM matrix should be symmetric, however when it's implemented 
     % for parallel computing, the imput data is just a subset of the original data.
     % so we didn't adjust the ssM matrix to make it symmetric.
     %ssM{tt} = max(ssM{tt},ssM{tt}'); 
     %ssM{tt} = tril(ssM{tt},-1) + triu(ssM{tt},1);
     %ssM{tt} = squareform(ssM{tt});
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
test = [ones(length(IND),1)*length(x) a(IND)+b(IND) a(IND)];
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


