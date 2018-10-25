function matsnpinfo(matlabFile,outputFile)

% FUNCTION snpinfo(matlabFile)

% SNPINFO gets SNP information from SNP data mat-file
%
% INPUTS:
%   matlabFile - a .mat file consists a structure array SNPdata with
%   the following fields:
%         rsid:snp names
%         data:genotype data
%         chr: chromosome id
%         loc: physical location
%         pheno: sample's phenotype
%         fid: sample's family id
%         pid: sample id
%         gender: sample gender
%
% OUTPUTs:
%   outputFile - a text file with 4 columns: SNP chromosome id, rs number, 
%                start location (all zeros), end location

load(matlabFile)

chr = SNPdata.chr;
rsid = SNPdata.rsid;
loc1 = zeros(length(chr),1);
loc2 = SNPdata.loc;

output = table(chr,rsid,loc1,loc2);
writetable(output,outputFile,'Delimiter',' ','WriteVariableNames',0);
