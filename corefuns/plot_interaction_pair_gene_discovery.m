function [ind1_gene_unique_new ind2_gene_unique_new heat_ssM_new] = plot_interaction_pair_gene_discovery(pathname1,pathname2,snppathwayfile,option,resultfile)
addpath /project/csbio/wwang/matlab_functions/misc/
addpath /project/csbio/wwang/matlab_functions/export_fig

%pathname1 = 'KEGG_RIBOSOME';
%pathname2 = 'KEGG_PARKINSONS_DISEASE';
%option = 2; % 1: protective; 2: risk
%resultfile = 'BPM_chi2_density0.1_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat';
%snppathwayfile='snp_pathway_min10_max300.mat';

load(resultfile)
ssM_dis = squareform(ssM{option});

clear chi2* ssM

load BPMind.mat
load SNPdata_mm.mat
load(snppathwayfile)
ii = find(ismember(pathwayNames,pathname1)==1);
jj = find(ismember(pathwayNames,pathname2)==1);
A = zeros(length(pathwayNames),length(pathwayNames));
A(ii,jj) = 1;
A(jj,ii) = 1;
A = squareform(A);
nn = find(A==1);

ind1 = BPMind1{nn};
ind2 = BPMind2{nn};

ssM_all_dis = ssM_dis;

ssM_dis = ssM_dis(ind1,ind2);
ind1_snp = snpRSNums(ind1);
ind2_snp = snpRSNums(ind2);

ind1 = 1:length(ind1);
ind2 = 1:length(ind2);

load snpgenemapping_50kb.mat
load /project/chadm/wwang/BPM2015/sourcedata/c2.all.v3.0.symbols_CP.mat

ind1_gp_idx = find(gp_matrix(:,find(ismember(pathwayNames,pathname1)==1))==1);
ind2_gp_idx = find(gp_matrix(:,find(ismember(pathwayNames,pathname2)==1))==1);

for i=1:length(ind1)
	ind1_gene_tmp = geneList(find(snp_gene_matrix(find(ismember(snpList,ind1_snp(i))==1),:)==1));
	ind1_gene_tmp = intersect(ind1_gene_tmp,uniquenames(ind1_gp_idx)'); 
	if length(ind1_gene_tmp)==1
		ind1_gene(i) = ind1_gene_tmp;
	elseif length(ind1_gene_tmp)>1
		ind1_gene(i) = cellstr(strjoin(ind1_gene_tmp','/'));
	else
		ind1_gene(i) = ind1_gene(i-1); % specical arrange
	end
end

for i=1:length(ind2)
	ind2_gene_tmp = geneList(find(snp_gene_matrix(find(ismember(snpList,ind2_snp(i))==1),:)==1));
	ind2_gene_tmp = intersect(ind2_gene_tmp,uniquenames(ind2_gp_idx)');
        if length(ind2_gene_tmp)==1
                ind2_gene(i) = ind2_gene_tmp;
        elseif length(ind2_gene_tmp)>1
                ind2_gene(i) = cellstr(strjoin(ind2_gene_tmp','/'));
	else
		ind2_gene(i) = cellstr('');
        end
end

heat_ssM_dis = zeros(length(unique(ind1_gene)),length(unique(ind2_gene)));

ind1_gene_unique = unique(ind1_gene);
ind2_gene_unique = unique(ind2_gene);

for i=1:length(ind1_gene_unique);
	ind1_idx{i} = find(ismember(ind1_gene,ind1_gene_unique{i})==1);
end

for i=1:length(ind2_gene_unique);
        ind2_idx{i} = find(ismember(ind2_gene,ind2_gene_unique{i})==1);
end

for i=1:length(ind1_gene_unique);
	for j=1:length(ind2_gene_unique)
		heat_ssM_dis(i,j) = nnz(ssM_dis(ind1(ind1_idx{i}),ind2(ind2_idx{j})));
	end
end
	
heat_ssM = heat_ssM_dis>0;
% plot entire heatmap

[tmp rowidx1] = sort(sum(heat_ssM,2),'descend');
[tmp colidx1] = sort(sum(heat_ssM,1),'descend');

heat_ssM_new = heat_ssM(rowidx1,colidx1);

ind1_gene_unique_new = ind1_gene_unique(rowidx1);
ind2_gene_unique_new = ind2_gene_unique(colidx1);

heatmap(heat_ssM_new,ind2_gene_unique_new,ind1_gene_unique_new,[],'TickAngle',45,'ShowAllTicks',1,'GridLines', '-','Colormap',[1,1,1; 1,0,0]);
setfig
export_fig('heatmap_all_rw.pdf')
close all

heatmap(heat_ssM_new,ind2_gene_unique_new,ind1_gene_unique_new,[],'TickAngle',45,'ShowAllTicks',1,'GridLines', '-','Colormap',[0,0,0; 1,0,0]);
setfig
export_fig('heatmap_all_rb.pdf')

close all

% rowidx1 = find(sum(heat_ssM,2)>(size(heat_ssM,2)*0.2));
% colidx1 = find(sum(heat_ssM,1)>(size(heat_ssM,1)*0.2));

rown1 = nnz(sum(heat_ssM_new,2)>(size(heat_ssM_new,2)*0.2));
coln1 = nnz(sum(heat_ssM_new,1)>(size(heat_ssM_new,1)*0.2));
%rowidx1 = find(sum(heat_ssM_new,2)>(size(heat_ssM_new,2)*0.2));
%colidx1 = find(sum(heat_ssM_new,1)>(size(heat_ssM_new,1)*0.2));
rowidx1 = 1:rown1;
colidx1 = 1:coln1;

heat_ssM_sub = heat_ssM_new(rowidx1,colidx1);

ind1_gene_unique_sub = ind1_gene_unique_new(rowidx1);
ind2_gene_unique_sub = ind2_gene_unique_new(colidx1);

% f = clustergram(heat_ssM_sub*1,'standardize',3);
% tmp = f.ColumnLabels;
% colidx = cellfun(@(x)str2num(x),tmp);

% tmp = f.RowLabels;
% rowidx = cellfun(@(x)str2num(x),tmp)';

figure
% heatmap(M(ind1,ind2),glist2u,glist1u,[],'TickAngle',45,'ShowAllTicks',true,'GridLines', '-');
heatmap(heat_ssM_new(rowidx1,colidx1),ind2_gene_unique_sub,ind1_gene_unique_sub,[],'TickAngle',45,'ShowAllTicks',1,'GridLines', '-','Colormap',[1,1,1; 1,0,0]);
setfig
export_fig('heatmap_rw.pdf')

close all

heatmap(heat_ssM(rowidx1,colidx1),ind2_gene_unique_sub,ind1_gene_unique_sub,[],'TickAngle',45,'ShowAllTicks',1,'GridLines', '-','Colormap',[0,0,0; 1,0,0]);
setfig
export_fig('heatmap_rb.pdf')







