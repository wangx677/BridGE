function geneset_of_interest(genes,ssmFile,snpgenemapping,outputname,samplePerm)

% this script is used to identify pathways and genes that are interacting with the given geneset
% genes: list of genes in a .mat file
% ssmFile: SNP-SNP interaction file
% snpgenemapping:
% outputname:
% binaryNetwork:
% snpPerms:
% minPath:
% netDiensity:

load(genes)
load(snpgenemapping)

% get snp index that can be mapped to genes
gene_ind = find(ismember(snp2gene.genelist,genes));
gene_ind = find(sum(snp2gene.sgmatrix(:,gene_ind),2)~=0);

load SNPdataAR.mat
targt_snps = SNPdata.rsid(gene_ind);

% check each SNP's interaction with the input genes to see if they're more than expected
% using ranksum test
non_gene_ind = setdiff(1:length(SNPdata.rsid),gene_ind);

for i = 0:samplePerm
     load(sprintf('%s_R%s.mat',ssmFile,num2str(i)))
     [p q] = size(ssM{1});
     if min(p,q) == 1;
          for tt=1:2
               ssM{tt} = squareform(ssM{tt});
          end
     end

     for tt=1:2
          for j = 1:length(SNPdata.rsid)
               snp_gi_enrich_ranksum{tt}(i+1,j) = ranksum(ssM{tt}(gene_ind,j),ssM{tt}(non_gene_ind,j),'tail','right');
          end

          a = length(SNPdata.rsid);
          b = length(gene_ind); 
          c = sum(ssM{tt}(gene_ind,:)>=0.2);
          d = sum(ssM{tt}>=0.2);
          
          snp_gi_enrich_hyge{tt}(i+1,:) = hygetest(a,b,c,d);
     end
     
end

for tt=1:2
     snp_gi_enrich_ranksum{tt} = -log10(snp_gi_enrich_ranksum{tt});
end

% get a set of candiate SNPs with signficance
 
outputFile = sprintf('geneset_of_interest_%s_%s.mat',outputname,ssmFile);
save(outputFile,'snp_gi_enrich*','genes','-v7.3')
