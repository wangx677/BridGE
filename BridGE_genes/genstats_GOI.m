function [bpm_local bpm_local_pv density_bpm density_bpm_expected wpm_local wpm_local_pv density_wpm density_wpm_expected denseidx path_degree path_degree_pv] = genstats(ssmFile,bpmindFile,binaryNetwork,snpPerms,minPath,netDensity)

% FUNCTION [bpm_local density_bpm density_bpm_expected wpm_local density_wpm density_wpm_expected denseidx path_degree] = genstats(ssmFile,bpmindFile,minPath,binaryNetwork,netDensity) 
%
% GENSTATS calculates general statistics for BPM (chi-square global,chi-square local, density), WPM (chi-square global,chi-square local, density) and 
%    PATH (ranksum)
%
% INPUTS:
%   ssmFile - SNP-SNP interaction file name
%   bpmindFile - a mat-file that includes SNP indexes of pathways for all possbile BPMs, WPMs
%   binaryNetwork - 1: binary network; 0: weighted network
%   snpPerms - number of SNP permutation
%   netDensity - network binarization density based on considering both protecitve/risk interactions
% 
% OUTPUTS:
%   An output file named genstats_<ssmFile>.mat or genstats_<ssmFile>_density<netDensity>.mat
%   

if (nargin > 6)
    error('Incorrect number of input arguments!');
end

load(bpmindFile)

load(sprintf('%s.mat',ssmFile))
[p q] = size(ssM{1});
if min(p,q) == 1;
     for tt=1:2
          ssM{tt} = squareform(ssM{tt});
     end
end

[p q] = size(ssM{1});
if binaryNetwork==1 & exist('netDensity')==1
     tmp = max(ssM{1},ssM{2});
     tmpcut = quantile(tmp(:),1-netDensity);
     for tt=1:2
          ssM{tt} = ssM{tt}>=tmpcut;
     end
     
     clear tmp tmpcut
elseif binaryNetwork==1 & exist('netDensity')~=1
     for tt=1:2
          ssM{tt} = ssM{tt}>0;
     end
     netDensity = nnz(max(ssM{1},ssM{2}))/(p*q)
end

% analyze protective and risk effect individually 
bpmsize = BPM.size;
wpmsize = WPM.size;
pathsize = WPM.indsize;
for tt=1:2
   [bpm_local{tt} bpm_local_pv{tt} density_bpm{tt} density_bpm_expected{tt} denseidx{tt} wpm_local{tt} wpm_local_pv{tt} density_wpm{tt} density_wpm_expected{tt} path_degree{tt} path_degree_pv{tt}] = rungenstats(full(ssM{tt}),BPM.ind1,BPM.ind2,BPM.ind1size,BPM.ind2size,bpmsize,binaryNetwork,snpPerms,minPath,WPM.ind,wpmsize,pathsize);
end

outputFile = sprintf('genstats_%s.mat',ssmFile);

if binaryNetwork==1
     save(outputFile,'bpm_local*','wpm_local*','path_degree*','denseidx','density_*','netDensity','-v7.3')
else
     save(outputFile,'bpm_local*','wpm_local*','path_degree*','denseidx','density_*','-v7.3')
end

end
