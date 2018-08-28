function data = bindataar(dataFile);

% FUNCTION data = bindataar(dataFile);
%
% BINDATAAR binarize 012 format SNP data based on recessive assumption.
%
% INPUTS:
%   dataFile - name of the data file. This .mat file consists 
%     a structure array SNPdata with the following fields:
%     - rsid:snp names
%     - data:genotype data
%     - chr: chromosome id 
%     - loc: physical location
%     - pheno: sample's phenotype
%     - fid: sample's family id
%     - pid: sample id
%     - gender
%
%  OUTPUTS:
%   A mat-file SNPdataAR.mat
%   

load(dataFile)
SNPdata.data(SNPdata.data==1) = 0;
SNPdata.data(SNPdata.data==2) = 1;
save('SNPdataAR.mat','SNPdata', '-v7.3')

data = SNPdata.data;

end
