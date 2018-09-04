% convert orf names to gene symbols in genesets

cd /project/csbio/wwang/BridGE/refdata

orf2gene = readtable('/project/csbio/wwang/BridGE/refdata/orf2gene.txt','FileType','text','ReadVariableNames',false);

file={'Costanzo2016_S6_HIER_pathway.mat','Costanzo2016_S6_HIER_process.mat','Costanzo2016_S12_complexes.mat','Costanzo2016_S6_HIER_compartment.mat','Costanzo2016_S5_SAFE.mat'};

for i=1:length(file)
     load(file{i})
     geneset.orfnames = geneset.genenames;
     for j=1:length(geneset.orfnames)
          ind = find(ismember(orf2gene.Var1,geneset.orfnames{j}));
          if isempty(ind)
               geneset.genenames(j) = geneset.orfnames(j);
          else
               geneset.genenames(j) = orf2gene.Var2(ind);
          end
     end
     
     save(file{i},'geneset');
     clear geneset
end

% convert orf names to gene symbols in snp to gene mapping file

orf2gene = readtable('/project/csbio/wwang/BridGE/refdata/orf2gene.txt','FileType','text','ReadVariableNames',false);

dirs = dir('*yeast*');
for i=1:length(dirs)
     cd(sprintf('/project/csbio/wwang/BridGE/%s',dirs(i).name))
     if  exist('snpgenemapping_500b.mat','file')==2
          system('mv snpgenemapping_500b.mat snpgenemapping_500bp.mat')
     end
     load('snpgenemapping_500bp.mat')
     snp2gene.orflist = snp2gene.genelist;
     for j=1:length(snp2gene.orflist)
          ind = find(ismember(orf2gene.Var1,snp2gene.orflist{j}));
          if isempty(ind)
               snp2gene.genelist(j) = snp2gene.orflist(j);
          else
               snp2gene.genelist(j) = orf2gene.Var2(ind);
          end
     end

     save('snpgenemapping_500bp.mat','snp2gene')
     clear snp2gene
end

% add geneset source 
dirs = dir('*yeast*GI*');
for i=1:length(dirs)
     cd(sprintf('/project/csbio/wwang/BridGE/%s',dirs(i).name))
     load('snp_pathway_min5_max300.mat')
     snpset.genesets = 'Costanzo2016_S6_HIER_pathway.mat';
     save('snp_pathway_min5_max300.mat','snpset');
end


     
