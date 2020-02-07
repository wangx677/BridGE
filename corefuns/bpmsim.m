function sim = bpmsim(BPMind1x,BPMind2x,BPMind1y,BPMind2y,option,ssM)

% FUNCTION sim = bpmsim(BPMind1x,BPMind2x,BPMind1y,BPMind2y,option,ssM)
%
% BPMSIM computes pairwise similarity for the given two sets of BPMs (X and Y).
% It also can be used to compute pairwise similarity for a given set (X=Y).
%
% INPUTS:
%    BPMind1x - BPMind1 for BPM set X
%    BPMind2x - BPMind2 for BPM set X
%    BPMind1y - BPMind1 for BPM set Y
%    BPMind2y - BPMind2 for BPM set Y
%    option - option=1 calculates overlap based on interactions;
%         option=2 calculates overlap based on SNPs
%    ssM - binarized SNP-SNP interaction matrix. It is required when option=1
% OUTPUTS:
%    sim - pairwise similarity matrix (overlap coefficient)

t1 = length(BPMind1x);
t2 = length(BPMind1y);

% check if BPMind1 is equal to BPMind2
TEST = isequal(BPMind1x,BPMind1y)==1 & isequal(BPMind2x,BPMind2y)==1;

BPMind1x = repmat(BPMind1x',1,t2);
BPMind2x = repmat(BPMind2x',1,t2);

BPMind1y = repmat(BPMind1y,t1,1);
BPMind2y = repmat(BPMind2y,t1,1);

if TEST==1
     BPMind1x = cellsquareform(BPMind1x);
     BPMind2x = cellsquareform(BPMind2x);
     BPMind1y = cellsquareform(BPMind1y);
     BPMind2y = cellsquareform(BPMind2y);
end

if option==1

	mm =cellfun(@(x,y,z,r) min(nnz(ssM(x,y)),nnz(ssM(z,r))),BPMind1x,BPMind2x,BPMind1y,BPMind2y);

	ind1 = cellfun(@(x,y) intersect(x,y),BPMind1x,BPMind1y,'UniformOutput',false);
	ind2 = cellfun(@(x,y) intersect(x,y),BPMind2x,BPMind2y,'UniformOutput',false);
	
	simtmp1 = cellfun(@(x,y,z) nnz(ssM(x,y)),ind1,ind2);
	simtmp1 = simtmp1./mm;

     ind1 = cellfun(@(x,y) intersect(x,y),BPMind1x,BPMind2y,'UniformOutput',false);
     ind2 = cellfun(@(x,y) intersect(x,y),BPMind2x,BPMind1y,'UniformOutput',false);

     simtmp2 = cellfun(@(x,y,z) nnz(ssM(x,y)),ind1,ind2);
	simtmp2 = simtmp2./mm;	

elseif option==2

	mm = cellfun(@(x,y,z,r) min(length(x)*length(y),length(z)*length(r)),BPMind1x,BPMind2x,BPMind1y,BPMind2y);

	ind1 = cellfun(@(x,y) intersect(x,y),BPMind1x,BPMind1y,'UniformOutput',false);
     ind2 = cellfun(@(x,y) intersect(x,y),BPMind2x,BPMind2y,'UniformOutput',false);

	simtmp1 = cellfun(@(x,y) (length(x)*length(y)),ind1,ind2);
	simtmp1 = simtmp1./mm;

	ind1 = cellfun(@(x,y) intersect(x,y),BPMind1x,BPMind2y,'UniformOutput',false);
     ind2 = cellfun(@(x,y) intersect(x,y),BPMind2x,BPMind1y,'UniformOutput',false);

	simtmp2 = cellfun(@(x,y) (length(x)*length(y)),ind1,ind2);
	simtmp2 = simtmp2./mm;
end

sim = max(simtmp1,simtmp2);

if TEST==1
     sim = squareform(sim);
end

