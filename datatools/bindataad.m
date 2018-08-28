function data = bindataad(dataFile);

% FUNCTION data = bindataad(dataFile);
%
% BINDATAAD binarize 012 format SNP data based on dominant assumption.
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
%    - gender
% 
% OUTPUTS:
%   a mat-file SNPdataAD.mat
%   

load(dataFile)
SNPdata.data(SNPdata.data==1) = 1;
SNPdata.data(SNPdata.data==2) = 1;
save('SNPdataAD.mat','SNPdata', '-v7.3')

data = SNPdata.data;

end
