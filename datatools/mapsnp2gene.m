function mapsnp2gene(snpAnnotation,geneAnnotation,mappingDistance,option,outputFile)

% FUNCTION mapsnp2gene(snpAnnotation,geneAnnotation,mappingDistance,option,outputFile)
%
% SNP2GENE  Map SNPs to genes based on given criteria
%
% INPUTS:
% snpAnnotation: plink bim file, or file in same format
% geneAnnotation: a text file in the following format: 
%  one row per gene, each row has chromosome, start and stop positions 
%  (base-pair) and then gene name. For example,
%  7 20140803 20223538 7A5
%  19 63549983 63556677 A1BG
%  10 52236330 52315441 A1CF
%  8 43266741 43337485 A26A1
% mappingDistance: mapping distance
% option: 'snplist','matrix'
%   'snplist': output is a list of SNPs in text file
%   'matrix': output is Mat-file with a structure array
%       including the following fields
%       snplist: list of SNPs
%       genelist: list of genes
%       sgmatrix: a binary matrix with rows are SNPs and columns are genes 
%
% OUTPUTs:
%   Depends on the value of input 'option'.
%

% read data from input files
[snpchr snprs tmp1 snploc tmp2 tmp3] = textread(snpAnnotation,'%s%s%s%d%s%s');
clear tmp1 tmp2 tmp3
[genechr geneloc1 geneloc2 genes] = textread(geneAnnotation,'%s%d%d%s');


% create empty arrays
if isequal(option,'snplist')
     snplist=[];
elseif isequal(option,'matrix')
     snprstmp = [];
     sgmatrix = [];
end

% loop from chromosome 1 to 22 to find mapping between SNPs and genes based
% on given distance
for i=1:22
     % find SNPs and genes that are located on same chromosome	
     indsnp = find(ismember(snpchr,num2str(i))==1);
     indgene = find(ismember(genechr,num2str(i))==1);

     % use repmat to create a snp location matrix. Rows are SNPs,
     % columns are repeated SNP locations and the number 
     % of columns is equal to the number of genes	
     snpmat = repmat(snploc(indsnp),1,length(indgene));

     % use repmat to create gene start/end location matrix with extanded 
     % window size to match the mapping distance. Columns are genes, 
     % rows are repeated gene start/end locations, and the number 
     % of rows is equal to the number of SNPs
     genemat1 = repmat(geneloc1(indgene)',length(indsnp),1) - mappingDistance;
     genemat2 = repmat(geneloc2(indgene)',length(indsnp),1) + mappingDistance;

     % create a sparse matrix for SNP and gene mapping
     snp2gene = sparse(length(indsnp),length(genes));

     % assign SNP-gene pairs relationship
     snp2gene(:,indgene) = sparse((snpmat >= genemat1) & (snpmat <= genemat2));

    % attach the loop results to general results
     if isequal(option,'snplist')
          ind = find(sum(snp2gene,2)~=0);
          snplist = [snplist;snprs(indsnp(ind))];
     elseif  isequal(option,'matrix')
          snprstmp = [snprstmp; snprs(indsnp)];
          sgmatrix = [sgmatrix;snp2gene];
     end
end
     
clear snp2gene
                              
if isequal(option,'snplist')
     % output is a list of SNPs that can be mapped to gene sets
     fwritecell(outputFile,'%s',snplist);
elseif isequal(option,'matrix')
     % output includes SNP list, gene list, and a binary matrix
     % indsnp = find(sum(sgmatrix,2) ~= 0);
     snplist = snprstmp;
     indgene = find(sum(sgmatrix,1) ~= 0);
     genelist = genes(indgene);
     sgmatrix = sgmatrix(:,indgene);
     snp2gene.snplist = snplist;
     snp2gene.genelist = genelist;
     snp2gene.sgmatrix = sgmatrix;

     % make sure the order of snp2gene.snplist is the same as the snpAnnotation file
     if isequal(snp2gene.snplist,snprs)~=1
          [tmp ida] = ismember(snprs,snp2gene.snplist);
          snp2gene.snplist = snp2gene.snplist(ida);
          snp2gene.sgmatrix = snp2gene.sgmatrix(ida,:);
     end
     save(outputFile,'snp2gene','-v7.3')
end


