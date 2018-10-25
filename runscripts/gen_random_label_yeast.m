groups = readtable('/project/csbio/wwang/yeast_1011strains/maf005_groups.txt','FileType','text','ReadVariableNames',false);
groups.Var2(find(cellfun(@(x)isempty(x),table2array(groups(:,2)))))={'other'};
uniquegroups = unique(groups.Var2);

rng(100)

load SNPdataAR.mat
rand_pheno = 100*ones(length(SNPdata.pheno),100);
for i=1:100
     for j=1:length(uniquegroups)
          ind = find(ismember(groups.Var2,uniquegroups{j}));
          ind1 = find(ismember(SNPdata.fid,groups.Var1(ind))); 
          if isempty(ind1)~=1
               ind2 = randperm(length(ind1));
               rand_pheno(ind1,i) = SNPdata.pheno(ind1(ind2));
          end
      end
end               

load SNPdataAR.mat
SNPdata.rand_pheno = rand_pheno;
save SNPdataAR.mat SNPdata

load SNPdataAD.mat
SNPdata.rand_pheno = rand_pheno;
save SNPdataAD.mat SNPdata

load gwas_data_final.mat
SNPdata.rand_pheno = rand_pheno;
save gwas_data_final.mat SNPdata
