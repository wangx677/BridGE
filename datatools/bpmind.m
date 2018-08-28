function bpmind(snpPathwayFile)

% FUNCTION bpmind(snpPathwayFile)
%
% BPMIND get the snp indexes of pathways for all possbile BPMs, WPMs
%
% INPUTS:
%   snpPathwayFile - a mat-file that includes a structure array snpset that
%     includes the following fields: snpPathwayMatrix, pathwaynames, snpList, and PathSize
%
% OUTPUTS:
%  BPM - a structure array that includes the following fields 
%  WPM - a structure array that includes the following fields
%

load(snpPathwayFile);

snp_pathway = snpset.spmatrix;
pp = size(snp_pathway,2);
BPMind1 = cell(pp,pp);
BPMind2 = cell(pp,pp);
WPMind = cell(1,pp);
path1idx = zeros(pp,pp);
path2idx = zeros(pp,pp);

for ii=1:pp
        for jj=ii+1:pp
                BPMind1{jj,ii} = find(snp_pathway(:,ii)==1 & snp_pathway(:,jj)~=1);
                BPMind2{jj,ii} = find(snp_pathway(:,ii)~=1 & snp_pathway(:,jj)==1);
		path1idx(jj,ii) = ii;
		path2idx(jj,ii) = jj;
        end
        WPMind{ii} = find(snp_pathway(:,ii)==1);
end

BPM.path1idx = squareform(path1idx);
BPM.path2idx = squareform(path2idx);
BPM.ind1 = cellsquareform(BPMind1);
BPM.ind2 = cellsquareform(BPMind2);
BPM.ind1size = cellfun(@(x) length(x),BPM.ind1);
BPM.ind2size = cellfun(@(x) length(x),BPM.ind2);
BPM.size = BPM.ind1size.*BPM.ind2size;

WPM.ind = WPMind;
WPM.indsize = cellfun(@(x) length(x),WPMind);
WPM.pathway = snpset.pathwaynames;
WPM.size = WPM.indsize.*WPM.indsize-WPM.indsize;

save('BPMind.mat','BPM','WPM','-v7.3');
