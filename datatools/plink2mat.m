function plink2mat(plinkRawFile,plinkBimFile,plinkFamFile,outputFile)

% FUNCTION plink2mat(plinkRawFile,plinkBimFile,plinkFamFile,outputFile)
% 
% PLINK2MAT convert plink .raw file to mat-file format
%
% This function converts .raw plink file into .mat file.
% It extracts all information from the .raw file and separates 
% the genotype information from the rest. It saves all into
% a <outputFile>.mat file. 
%
% INPUTS: 
% plinkRawFile - plink .raw file
% plinkBimFile - plink.bim file that associated with plinkRawFile
% plinkFamFile - plink.fam file that associated with plinkRawFile
% outputFile - name for output mat-file
%
% OUTPUTS:
% <outputFile>.mat
%   The .mat file consists a structure array SNPdata with
%   the following fields:
%   - rsid:snp names
%   - data:genotype data
%   - chr: chromosome id 
%   - loc: physical location
%   - pheno: sample's phenotype
%   - fid: sample's family id
%   - pid: sample id
%   - gender: sample gender


% read plink raw data
data = readcsv(plinkRawFile,' ');

tmpfid = data(:,1);
tmpfid = tmpfid(2:end);

tmppid = data(:,2);
tmppid = tmppid(2:end);

tmppheno = data(:,6);
try
     tmppheno = str2num(cell2mat(tmppheno(2:end)));
catch
     tmppheno = cellfun(@(x)str2num(x),tmppheno(2:end));
end

tmpgender = data(:,5);
tmpgender = str2num(cell2mat(tmpgender(2:end)));

data(:,1:6) = [];

tmprsid = data(1,:)';
A = regexp (tmprsid, '_', 'split');
tmprsid = cellfun(@(x) x{1},A,'uniform',0);

data(1,:) = [];

data = strrep(data,'NA','-1');

for i=1:size(data,1)
     tmpdata = data(i,:)';
     tmp = sprintf('%s*',tmpdata{:});
     trans_data(i,:) = sscanf(tmp, '%d*');
end

data = trans_data; 
clear trans_data

[chr rsid tmp1 loc tmp2 tmp3] = textread(plinkBimFile,'%d%s%s%d%s%s');
[fid pid tmp1 tmp2 gender pheno] = textread(plinkFamFile,'%s%s%s%s%d%f');

if isequal(tmprsid,rsid)~=1
     sprintf('rsid of plinkBimFile and matlab file does not match, please check!')
     exit
end

if isequal(tmpfid,fid)~=1 | isequal(tmppid,pid)~=1 | isequal(tmppheno,pheno)~=1 | isequal(tmpgender,gender)~=1
     sprintf('subjects of plinkBimFile and matlab file does not match, please check!')
     exit
end     

data = imputesnp(data);

if (isequal(unique(pheno),[1 2]'))
     pheno = pheno -1;
end

SNPdata.data = data;
SNPdata.rsid = rsid;
SNPdata.chr = chr;
SNPdata.loc = loc;
SNPdata.pheno = pheno;
SNPdata.fid = fid;
SNPdata.pid = pid;
SNPdata.gender = gender;

save(outputFile,'SNPdata','-v7.3')
