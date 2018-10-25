function check_gi_case_control_distribution(file,model,GIcutoff,outputfile)

% file = '/project/csbio/wwang/BridGE/project_yeast_SC4NQO01ugml_38h_complex_t25_b25_mhygeSSI/BPM_WPM_info/BPM8_DD_protective_fdr0.00_pv2.00e-04_rs10.13.xls';
% GIcutoff = 2
% outputfile

% file - BPM or WPM info file

% get yeast strain's group information
groups = readtable('/project/csbio/wwang/yeast_1011strains/maf005_groups.txt','FileType','text','ReadVariableNames',false);
groups.Var2(find(ismember(groups.Var2,''))) = {'other'};

% read GI pairs that support the BPM/WPM interaction
data = readtable(file,'Sheet',1);

% only evaluate GI paris that pass the GIcutoff threshold
ind = find(data.GI>GIcutoff);

% read genetic data
if strcmp(model,'DD')
     load('../SNPdataAD.mat');
elseif strcmp(model,'RR')
      load('../SNPdataAR.mat');
end

ncase = nnz(SNPdata.pheno==1);
ncontrol = nnz(SNPdata.pheno==0);
for i = 1:length(ind)
     % identify SNP pair
     snp1 = data.snps1(ind(i));
     snp2 = data.snps2(ind(i));
     
     % locate SNP pair genotype information
     ind1 = find(ismember(SNPdata.rsid,snp1));
     ind2 = find(ismember(SNPdata.rsid,snp2));
     freq00_control(i) = nnz(SNPdata.data(:,ind1)==0 & SNPdata.data(:,ind2)==0 & SNPdata.pheno==0)/ncontrol;
     freq00_case(i) = nnz(SNPdata.data(:,ind1)==0 & SNPdata.data(:,ind2)==0 & SNPdata.pheno==1)/ncase;
          
     freq01_control(i) = nnz(SNPdata.data(:,ind1)==0 & SNPdata.data(:,ind2)==1 & SNPdata.pheno==0)/ncontrol;
     freq01_case(i) = nnz(SNPdata.data(:,ind1)==0 & SNPdata.data(:,ind2)==1 & SNPdata.pheno==1)/ncase;

     freq10_control(i) = nnz(SNPdata.data(:,ind1)==1 & SNPdata.data(:,ind2)==0 & SNPdata.pheno==0)/ncontrol;
     freq10_case(i) = nnz(SNPdata.data(:,ind1)==1 & SNPdata.data(:,ind2)==0 & SNPdata.pheno==1)/ncase;

     freq11_control(i) = nnz(SNPdata.data(:,ind1)==1 & SNPdata.data(:,ind2)==1 & SNPdata.pheno==0)/ncontrol;
     freq11_case(i) = nnz(SNPdata.data(:,ind1)==1 & SNPdata.data(:,ind2)==1 & SNPdata.pheno==1)/ncase;

     % identify which group (case or control) has more 11 combinations      
     % group has higher 11 combinations is the group to evluate population    
     idx = find(SNPdata.data(:,ind1)==1 & SNPdata.data(:,ind2)==1);
     idx1 = find(SNPdata.data(:,ind1)==1 & SNPdata.data(:,ind2)==1 & SNPdata.pheno==0);
     idx2 = find(SNPdata.data(:,ind1)==1 & SNPdata.data(:,ind2)==1 & SNPdata.pheno==1);
 
     if length(idx2)>length(idx1)
          tmp = groups.Var2(ismember(groups.Var1,SNPdata.fid(idx2)));
     else
          tmp = groups.Var2(ismember(groups.Var1,SNPdata.fid(idx1)));
     end

     ugroup =  unique(tmp);
     ugroup_no = cellfun(@(x)nnz(ismember(tmp,x)),ugroup);
     [tmp idx] = sort(ugroup_no,'descend');
     ugroup = ugroup(idx);
     ugroup_no = ugroup_no(idx);
          
     for j=1:length(ugroup)
          tmp_ugroup = strsplit(ugroup{j},'.');
          tmp_ugroup = tmp_ugroup{1};
          tmp_percent = ugroup_no(j)/nnz(ismember(groups.Var2(find(ismember(groups.Var1,SNPdata.fid))),ugroup{j}));
          tmp_groupinfo{j} =  sprintf('cluster%s_%d_%.2f',tmp_ugroup,ugroup_no(j),tmp_percent); 
     end
     groupinfo{i} = strjoin(tmp_groupinfo,',');
     clear tmp*
end

output = [freq00_control' freq00_case' freq01_control' freq01_case' freq10_control' freq10_case' freq11_control' freq11_case'];
bar(output')
setfig
saveas(gcf,sprintf('%s.pdf',outputfile))
saveas(gcf,sprintf('%s.png',outputfile))
