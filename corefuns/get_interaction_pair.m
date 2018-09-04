function [output_path1_snp, output_path2_snp] = get_interaction_pair(pathname1,pathname2,effect,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,outputfile);

%pathname1 = 'KEGG_RIBOSOME';
%pathname2 = 'KEGG_PARKINSONS_DISEASE';
% effect = 'proective'; % 1: protective; 2: risk
%resultfile = 'BPM_chi2_density0.1_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat';
%snppathwayfile='snp_pathway_min10_max300.mat';

load(ssmfile)
load(bpmindfile)
load(snppathwayfile)
genesetfile = snpset.genesets;
load(sprintf('%s/refdata/%s',getenv('BRIDGEPATH'),genesetfile));
load(snpgenemappingfile)
load('SNPdataAR.mat'); SNPdataAR = SNPdata.data;
load('SNPdataAD.mat'); SNPdataAD = SNPdata.data;

model = strsplit(ssmfile,'_');
model = model{end-1};

if ismember(effect,'protective')
     ssM_dis = squareform(ssM{1});
     if exist('maxidx')
          maxidx = squareform(maxidx{1});
     end
else
     ssM_dis = squareform(ssM{2});
     if exist('maxidx')
          maxidx = squareform(maxidx{2});
     end
end
ssM_all_dis = ssM_dis;
N = size(ssM_dis,1);
clear ssM

if isequal(pathname1,pathname2)
     ii = find(ismember(snpset.pathwaynames,pathname1)==1);
     ind1 = WPM.ind{ii};
     ind2 = WPM.ind{ii};
     ssM_dis = ssM_dis(ind1,ind2);
     ssM_dis = tril(ssM_dis);
     if exist('maxidx','var')
          maxidx = maxidx(ind1,ind2);
          maxidx  = tril(maxidx);
     end
else
     ii = find(ismember(snpset.pathwaynames,pathname1)==1);
     jj = find(ismember(snpset.pathwaynames,pathname2)==1);
     A = zeros(length(snpset.pathwaynames),length(snpset.pathwaynames));
     A(ii,jj) = 1;
     A(jj,ii) = 1;
     A = squareform(A);
     nn = find(A==1);

     ind1 = BPM.ind1{nn};
     ind2 = BPM.ind2{nn};
     ssM_dis = ssM_dis(ind1,ind2);
     if exist('maxidx','var')
          maxidx = maxidx(ind1,ind2);
     end
end

ind1_snp = SNPdata.rsid(ind1);
ind2_snp = SNPdata.rsid(ind2);

SNPdataAR1 = SNPdataAR(:,ind1);
SNPdataAR2 = SNPdataAR(:,ind2);

SNPdataAD1 = SNPdataAD(:,ind1);
SNPdataAD2 = SNPdataAD(:,ind2);

ind1 = 1:length(ind1);
ind2 = 1:length(ind2);

ind1_gp_idx = find(geneset.gpmatrix(:,find(ismember(geneset.pathwaynames,pathname1)==1))==1);
ind2_gp_idx = find(geneset.gpmatrix(:,find(ismember(geneset.pathwaynames,pathname2)==1))==1);

for i=1:length(ind1)
	ind1_gene_tmp = snp2gene.genelist(find(snp2gene.sgmatrix(find(ismember(snp2gene.snplist,ind1_snp(i))==1),:)==1));
	ind1_gene_tmp = intersect(ind1_gene_tmp,geneset.genenames(ind1_gp_idx)'); 
	if length(ind1_gene_tmp)==1
		ind1_gene(i) = ind1_gene_tmp;
	elseif length(ind1_gene_tmp)>1
		ind1_gene(i) = cellstr(strjoin(ind1_gene_tmp','/'));
	else
		ind1_gene(i) = cellstr(''); % specical arrange
	end
end

if isequal(pathname1,pathname2)
     ind2_gene = ind1_gene;
else
     for i=1:length(ind2)
     	ind2_gene_tmp = snp2gene.genelist(find(snp2gene.sgmatrix(find(ismember(snp2gene.snplist,ind2_snp(i))==1),:)==1));
     	ind2_gene_tmp = intersect(ind2_gene_tmp,geneset.genenames(ind2_gp_idx)');
             if length(ind2_gene_tmp)==1
                     ind2_gene(i) = ind2_gene_tmp;
           elseif length(ind2_gene_tmp)>1
                     ind2_gene(i) = cellstr(strjoin(ind2_gene_tmp','/'));
	     else
	     	ind2_gene(i) = cellstr('');
           end
     end
end

[i j] = find(ssM_dis>0.2);

snps1 = ind1_snp(i);
snps2 = ind2_snp(j);
genes1 = ind1_gene(i)';
genes2 = ind2_gene(j)';

GI = zeros(length(i),1);
GItype = cell(length(i),1);
freq_case = zeros(length(i),1);
freq_control = zeros(length(i),1);

for k=1:length(i)
     GI(k) = ssM_dis(i(k),j(k));
     if strcmp(model,'combined')
          if maxidx(i(k),j(k))==1
               GItype{k} = 'recessive';
               freq_case(k) = nnz(SNPdataAR1(SNPdata.pheno==1,i(k)).*SNPdataAR2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
               freq_control(k) = nnz(SNPdataAR1(SNPdata.pheno==0,i(k)).*SNPdataAR2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);           
               snp_pair(:,k) = SNPdataAR1(:,i(k)).*SNPdataAR2(:,j(k));
          elseif  maxidx(i(k),j(k))==2
               GItype{k} = 'dominant';
               freq_case(k) = nnz(SNPdataAD1(SNPdata.pheno==1,i(k)).*SNPdataAD2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
               freq_control(k) = nnz(SNPdataAD1(SNPdata.pheno==0,i(k)).*SNPdataAD2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);
               snp_pair(:,k) = SNPdataAD1(:,i(k)).*SNPdataAD2(:,j(k));
          elseif maxidx(i(k),j(k))==3
               GItype{k} = 'recessive_dominant';
               freq_case1 = nnz(SNPdataAD1(SNPdata.pheno==1,i(k)).*SNPdataAR2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
               freq_control1 = nnz(SNPdataAD1(SNPdata.pheno==0,i(k)).*SNPdataAR2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);
               freq_case2 = nnz(SNPdataAR1(SNPdata.pheno==1,i(k)).*SNPdataAD2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
               freq_control2 = nnz(SNPdataAR1(SNPdata.pheno==0,i(k)).*SNPdataAD2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);
               freqR1 = freq_case1/freq_control1;
               freqR2 = freq_case2/freq_control2;     

               if ismember(effect,'protective')
                    if freqR1>freqR2
                         freq_case(k) = freq_case2;
                         freq_control(k) = freq_control2; 
                         snp_pair(:,k) = SNPdataAR1(:,i(k)).*SNPdataAD2(:,j(k));
                    else
                         freq_case(k) = freq_case1;
                         freq_control(k) = freq_control1;
                         snp_pair(:,k) = SNPdataAD1(:,i(k)).*SNPdataAR2(:,j(k));
                    end
               elseif ismember(effect,'risk')
                     if freqR1>freqR2
                         freq_case(k) = freq_case1;
                         freq_control(k) = freq_control1;
                         snp_pair(:,k) = SNPdataAD1(:,i(k)).*SNPdataAR2(:,j(k));
                    else
                         freq_case(k) = freq_case2;
                         freq_control(k) = freq_control2;
                         snp_pair(:,k) = SNPdataAR1(:,i(k)).*SNPdataAD2(:,j(k));
                    end
               end          

          end
     elseif strcmp(model,'RR')

          GItype{k} = 'recessive';
          freq_case(k) = nnz(SNPdataAR1(SNPdata.pheno==1,i(k)).*SNPdataAR2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
          freq_control(k) = nnz(SNPdataAR1(SNPdata.pheno==0,i(k)).*SNPdataAR2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);
          snp_pair(:,k) = SNPdataAR1(:,i(k)).*SNPdataAR2(:,j(k));

     elseif strcmp(model,'DD')

          GItype{k} = 'dominant';
          freq_case(k) = nnz(SNPdataAD1(SNPdata.pheno==1,i(k)).*SNPdataAD2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
          freq_control(k) = nnz(SNPdataAD1(SNPdata.pheno==0,i(k)).*SNPdataAD2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);
          snp_pair(:,k) = SNPdataAD1(:,i(k)).*SNPdataAD2(:,j(k));

     elseif strcmp(model,'RD')

           GItype{k} = 'recessive_dominant';
           freq_case1 = nnz(SNPdataAD1(SNPdata.pheno==1,i(k)).*SNPdataAR2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
           freq_control1 = nnz(SNPdataAD1(SNPdata.pheno==0,i(k)).*SNPdataAR2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);
           freq_case2 = nnz(SNPdataAR1(SNPdata.pheno==1,i(k)).*SNPdataAD2(SNPdata.pheno==1,j(k)))/nnz(SNPdata.pheno==1);
           freq_control2 = nnz(SNPdataAR1(SNPdata.pheno==0,i(k)).*SNPdataAD2(SNPdata.pheno==0,j(k)))/nnz(SNPdata.pheno==0);
           freqR1 = freq_case1/freq_control1;
           freqR2 = freq_case2/freq_control2;

           if ismember(effect,'protective')
                if freqR1>freqR2
                     freq_case(k) = freq_case2;
                     freq_control(k) = freq_control2;
                     snp_pair(:,k) = SNPdataAR1(:,i(k)).*SNPdataAD2(:,j(k));
                else
                     freq_case(k) = freq_case1;
                     freq_control(k) = freq_control1;
                     snp_pair(:,k) = SNPdataAD1(:,i(k)).*SNPdataAR2(:,j(k));
                end
           elseif ismember(effect,'risk')
                if freqR1>freqR2
                     freq_case(k) = freq_case1;
                     freq_control(k) = freq_control1;
                     snp_pair(:,k) = SNPdataAD1(:,i(k)).*SNPdataAR2(:,j(k));
                else
                     freq_case(k) = freq_case2;
                     freq_control(k) = freq_control2;
                     snp_pair(:,k) = SNPdataAR1(:,i(k)).*SNPdataAD2(:,j(k));
                end
           end
     end

end

output_pair = table(snps1,genes1,snps2,genes2,GItype,freq_case,freq_control,GI);
output_pair = sortrows(output_pair,{'GI'},{'descend'});

if exist('outputfile','var')
     writetable(output_pair,sprintf('%s.xls',outputfile),'Sheet',1);
end

ind1_mean_GI  = sum(ssM_dis>0.2,2)/length(ind2);

for i=1:length(ind1)
     idx(i) = find(ismember(SNPdata.rsid,ind1_snp(i)));
end
ind1_mean_GI_bg = sum(ssM_all_dis(idx,:)>0.2,2)/size(ssM_all_dis,1);

snps = ind1_snp;
genes = ind1_gene';
snp_mean_gi = ind1_mean_GI;
snp_mean_gi_bg = ind1_mean_GI_bg;
gi_fold = snp_mean_gi./snp_mean_gi_bg;
gi_hyge = hygetest(N,length(ind2),sum(ssM_dis>0.2,2),sum(ssM_all_dis(idx,:)>0.2,2));
gi_fold(isnan(gi_fold)) = 0;
output_path1_snp = table(snps,genes,snp_mean_gi,snp_mean_gi_bg,gi_fold,gi_hyge);

[tmp ind] = sort(gi_fold,'descend');
output_path1_snp = output_path1_snp(ind,:);
output_path1_snp = output_path1_snp(find(output_path1_snp.gi_fold>1.2),:);
output_path1_snp = sortrows(output_path1_snp,{'gi_hyge'},{'descend'});

if exist('outputfile','var')
     writetable(output_path1_snp,sprintf('%s.xls',outputfile),'Sheet',2);
end

clear idx

if isequal(pathname1,pathname2)~=1     
     ind2_mean_GI  = sum(ssM_dis>0.2)/length(ind1)';
     
     for i=1:length(ind2)
          idx(i) = find(ismember(SNPdata.rsid,ind2_snp(i)));
     end
     ind2_mean_GI_bg = sum(ssM_all_dis(:,idx)>0.2)/size(ssM_all_dis,1);

     snps = ind2_snp;
     genes = ind2_gene';
     snp_mean_gi = ind2_mean_GI;
     snp_mean_gi_bg = ind2_mean_GI_bg;
     gi_fold = snp_mean_gi./snp_mean_gi_bg;
     gi_hyge = hygetest(N,length(ind1), sum(ssM_dis>0.2),sum(ssM_all_dis(:,idx)>0.2));
     snp_mean_gi = snp_mean_gi';
     snp_mean_gi_bg = snp_mean_gi_bg';
     gi_fold = gi_fold';
     gi_hyge = gi_hyge';
     gi_fold(isnan(gi_fold)) = 0;
     output_path2_snp = table(snps,genes,snp_mean_gi,snp_mean_gi_bg,gi_fold,gi_hyge);

     [tmp ind] = sort(gi_fold,'descend');
     output_path2_snp = output_path2_snp(ind,:);
     output_path2_snp = output_path2_snp(find(output_path2_snp.gi_fold>1.2),:);
     output_path2_snp = sortrows(output_path2_snp,{'gi_hyge'},{'descend'});
      
     if exist('outputfile','var')
          writetable(output_path2_snp,sprintf('%s.xls',outputfile),'Sheet',3);
          interaction = 'BPM';
          output_interaction = table(interaction,pathname1,pathname2);
          writetable(output_interaction,sprintf('%s.xls',outputfile),'Sheet',4);
     end
     clear idx
else
     output_path2_snp = nan;
     if exist('outputfile','var')
          interaction = 'WPM';
          output_interaction = table(interaction,pathname1);
          writetable(output_interaction,sprintf('%s.xls',outputfile),'Sheet',3);
     end
end

if exist('outputfile','var')
     save(sprintf('%s.mat',outputfile),'output*','snp_pair')
end

