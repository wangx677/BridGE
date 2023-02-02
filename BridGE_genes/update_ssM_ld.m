function update_ssM_ld(ssMfile)

% this function is to remove interactions that based on LD structure
% 

% load SNP-SNP interaction file
load(ssMfile);

% exit program if ssM_ld exist
if exist('ssM_ld','var')==1
     exit;
end

% load gwas data
load('gwas_data_final.mat')

% load LD file
ld = readtable('plink.ld','filetype','text');
n = length(SNPdata.rsid);
ldtmp = zeros(n,n);

for i=1:size(ld,1)
     idx1 = find(ismember(SNPdata.rsid,ld.SNP_A{i}));
     idx2 = find(ismember(SNPdata.rsid,ld.SNP_B{i}));
     ldtmp(idx1,idx2) = ld.R2(i);
end

ldtmp = max(ldtmp,ldtmp');

% convert SSM vector to matrix
for tt=1:2
     ssM{tt} = squareform(ssM{tt});
end

% compute SNP-SNP similarity based on interaction
for tt=1:2
     ssMr2{tt} = squareform(1-pdist((ssM{tt}>=0.2)*1,'jaccard'));
end

% get interaction degree
% for SNPs in LD, we keep SNPs with larger  interaction degrees
for tt=1:2
     ssM_ld{tt} = ssM{tt};
     ssMdegree{tt} = sum(ssM{tt}>=0.2);
end

% set genetic interaction score to 0 for redundant SNPs in LD
uchr = unique(SNPdata.chr);
for tt=1:2
     kk=1
     while nnz(ssMr2{tt}>=0.5 & ldtmp>=0.2)>0
          % start from SNPs with highest interaction degree
          [tmp1 tmp2] = find(ssMr2{tt}>=0.5 & ldtmp>=0.2);
          ind = unique([tmp1;tmp2]);
          [tmp idx] = sort(ssMdegree{tt}(ind),'descend');
          for i=1:length(ind)-1
               if sum(ssM{tt}(ind(idx(i)),:)) ~= 0 % skip SNP if it doesn't have interactions
                    for j=i+1:length(ind)
                         % if two SNPs are in LD and their interaction profiles have >=0.5 similarity, 
                         % set the 2nd SNP interaction profile to be all 0s 
                         if (ldtmp(ind(idx(i)),ind(idx(j)))>=0.2  & ssMr2{tt}(ind(idx(i)),ind(idx(j)))>=0.5)
                              ssM_ld{tt}(ind(idx(j)),:) = 0;
                              ssM_ld{tt}(:,ind(idx(j))) = 0;
                         end
                    end
               end
          end
          ssMr2{tt} = squareform(1-pdist((ssM_ld{tt}>=0.2)*1,'jaccard'));
          ssM{tt} = ssM_ld{tt};
          ssMdegree{tt} = sum(ssM{tt}>=0.2);
          kk=kk+1
     end
end

for tt=1:2
     ssM{tt} = squareform(ssM_ld{tt});
end

tmp = strsplit(ssMfile,'ssM_');
outputfile = sprintf('ssM_ld_%s',tmp{2});

save(outputfile,'ssM')
