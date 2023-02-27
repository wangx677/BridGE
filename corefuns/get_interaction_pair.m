function [output_path1_snp, output_path2_snp] = get_interaction_pair(pathname1,pathname2,effect,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,densitycutoff,outputfile);

% Inputs:
% pathname1: first pathway name 
% pathname2: second pathway name, if pathname1 and pathbname2 are the same, it is a WPM
% effect: 'proective' or 'risk.  1: protective; 2: risk
% ssmfile: SNP-SNP interaction file
% bpmindfile: file that stores BPM and WPM indexes
% snppathwayfile: SNP to pathway mapping file
% snpgenemappingfile: SNP to gene mapping file
% densitycutoff: use the defined density to binarize the SNP-SNP interaction network and find driver SNPs and genes
% outputfile: name of the output file

% load files
load(ssmfile)
load(bpmindfile)
load(snppathwayfile)
genesetfile = snpset.genesets;
load(sprintf('%s/refdata/%s',getenv('BRIDGEPATH'),genesetfile));
load(snpgenemappingfile)
load('SNPdataAR.mat'); SNPdataAR = SNPdata.data;
load('SNPdataAD.mat'); SNPdataAD = SNPdata.data;

% identify disease model
model = strsplit(ssmfile,'_');
model = model{end-1};

% if disease model is not RR, DD, RD and combined, we don't know which disease model it is
if strcmp(model,'RR')~=1 && strcmp(model,'DD')~=1 && strcmp(model,'RD')~=1 && strcmp(model,'combined')~=1
   model = '';
end

% identify effect
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

% binarize the SNP-SNP interaction matrix
% ssM_dis_bin is the binarized matrix
if exist('densitycutoff','var') == 1
   dcutoff = quantile(squareform(ssM_dis),1-densitycutoff);
   ssM_dis_bin = ssM_dis > dcutoff;
else
    ssM_dis_bin = ssM_dis > 0.2;
end

ssM_all_dis = ssM_dis;  % SNP-SNP interaction matrix
ssM_all_dis_bin = ssM_dis_bin;

N = size(ssM_dis,1); % number of SNPs
clear ssM

if isequal(pathname1,pathname2)
     % for WPM, identify SNP-SNP interaction matrix
     ii = find(ismember(snpset.pathwaynames,pathname1)==1);
     ind1 = WPM.ind{ii};
     ind2 = WPM.ind{ii};
     ssM_dis = ssM_dis(ind1,ind2);
     ssM_dis_bin = ssM_dis_bin(ind1,ind2);
     % WPM is a symmetric matrix, only keep lower half interactions
     ssM_dis = tril(ssM_dis); 
     ssM_dis_bin = tril(ssM_dis_bin);
     if exist('maxidx','var')
          maxidx = maxidx(ind1,ind2);
          maxidx  = tril(maxidx);
     end
else
     % for BPM, identify SNP-SNP interaction matrix
     if isfield(BPM,'path1')
           nn = find(ismember(BPM.path1,pathname1)==1 & ismember(BPM.path2,pathname2)==1);
     else
          nn = find(ismember(WPM.pathway(BPM.path1idx),pathname1)==1 & ismember(WPM.pathway(BPM.path2idx),pathname2)==1);
     end

     ind1 = BPM.ind1{nn};
     ind2 = BPM.ind2{nn};
     ssM_dis = ssM_dis(ind1,ind2);
     ssM_dis_bin = ssM_dis_bin(ind1,ind2);
     if exist('maxidx','var')
          maxidx = maxidx(ind1,ind2);
     end
end

% identify corresponding genes
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


% identify SNP-SNP interaction pairs
[i j] = find(ssM_dis_bin>0);

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

if isempty(model)~=1
   output_pair = table(snps1,genes1,snps2,genes2,GItype,freq_case,freq_control,GI);
else
   output_pair = table(snps1,genes1,snps2,genes2,GI);
end
output_pair = sortrows(output_pair,{'GI'},{'descend'});

if exist('outputfile','var')
     writetable(output_pair,sprintf('%s.xls',outputfile),'Sheet',1);
end

% identify driver SNPs/genes
% make WPM SNP-SNP interaction matrix symmetric again
if isequal(pathname1,pathname2)==1
     ssM_dis_bin = max(ssM_dis_bin,ssM_dis_bin');
end

ind1_mean_GI  = sum(ssM_dis_bin>0, 2)/length(ind2);

for i=1:length(ind1)
     idx(i) = find(ismember(SNPdata.rsid,ind1_snp(i)));
end
ind1_mean_GI_bg = sum(ssM_all_dis_bin(idx,:)>0, 2)/size(ssM_all_dis_bin,1);

snps = ind1_snp;
genes = ind1_gene';
snp_mean_gi = ind1_mean_GI;
snp_mean_gi_bg = ind1_mean_GI_bg;
gi_fold = snp_mean_gi./snp_mean_gi_bg;
gi_hyge = hygetest(N,length(ind2),sum(ssM_dis_bin>0, 2),sum(ssM_all_dis_bin(idx,:)>0, 2));
gi_fold(isnan(gi_fold)) = 0;
output_path1_snp = table(snps,genes,snp_mean_gi,snp_mean_gi_bg,gi_fold,gi_hyge);

[tmp ind] = sort(gi_fold,'descend');
output_path1_snp = output_path1_snp(ind,:);
output_path1_snp = output_path1_snp(find(output_path1_snp.gi_fold>1.2 & output_path1_snp.gi_hyge>=-log10(0.05)),:);
output_path1_snp = sortrows(output_path1_snp,{'gi_fold'},{'descend'});

if exist('outputfile','var')
     writetable(output_path1_snp,sprintf('%s.xls',outputfile),'Sheet',2);
end

clear idx

if isequal(pathname1,pathname2)~=1     
     ind2_mean_GI  = sum(ssM_dis_bin>0)/length(ind1)';
     
     for i=1:length(ind2)
          idx(i) = find(ismember(SNPdata.rsid,ind2_snp(i)));
     end
     ind2_mean_GI_bg = sum(ssM_all_dis_bin(:,idx)>0)/size(ssM_all_dis_bin,1);

     snps = ind2_snp;
     genes = ind2_gene';
     snp_mean_gi = ind2_mean_GI;
     snp_mean_gi_bg = ind2_mean_GI_bg;
     gi_fold = snp_mean_gi./snp_mean_gi_bg;
     gi_hyge = hygetest(N,length(ind1), sum(ssM_dis_bin>0),sum(ssM_all_dis_bin(:,idx)>0));
     snp_mean_gi = snp_mean_gi';
     snp_mean_gi_bg = snp_mean_gi_bg';
     gi_fold = gi_fold';
     gi_hyge = gi_hyge';
     gi_fold(isnan(gi_fold)) = 0;
     output_path2_snp = table(snps,genes,snp_mean_gi,snp_mean_gi_bg,gi_fold,gi_hyge);

     [tmp ind] = sort(gi_fold,'descend');
     output_path2_snp = output_path2_snp(ind,:);
     output_path2_snp = output_path2_snp(find(output_path2_snp.gi_fold>1.2 & output_path2_snp.gi_hyge>=-log10(0.05)),:);
     output_path2_snp = sortrows(output_path2_snp,{'gi_fold'},{'descend'});
      
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

