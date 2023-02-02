load stats_single_association_and_pairwise_interaction.mat
load SNPdataAD.mat
SNPdata_AD = SNPdata;
N = length(SNPdata_AD.pheno);
AD_SNP_protective(1,:) = hygetest(N,nnz(SNPdata_AD.pheno==0),sum(SNPdata_AD.data(find(SNPdata_AD.pheno==0),:)),sum(SNPdata_AD.data));
AD_SNP_risk(1,:) = hygetest(N,nnz(SNPdata_AD.pheno==1),sum(SNPdata_AD.data(find(SNPdata_AD.pheno==1),:)),sum(SNPdata_AD.data));

load ssM_mhygeSSI_alpha10.05_alpha20.05_DD_R0.mat
for tt=1:2
     ssM{tt} = squareform(ssM{tt});
end

% single SNP association vs GI
scatter(AD_SNP_protective,sum(ssM{1}>0.5)'/length(AD_SNP_protective_pv(1,:)));
xlabel('Single SNP association -log10(p-value)')
ylabel('GI degree (fraction)')
saveas(gcf,'scatter_single_SNP_association_GI_degreee_mhygeSSI.png')

% single SNP association vs BPM
load genstats_ssM_mhygeSSI_alpha10.05_alpha20.05_DD_R0.mat
ind = find(bpm_local_pv{1}<=0.05);
load('BPMind.mat');
BPMind1 = BPM.ind1(ind);
BPMind2 = BPM.ind2(ind);

for i=1:length(AD_SNP_protective_pv(1,:))
     ind1 = find(cell2mat(cellfun(@(x)sum(ismember(x,i)>0),BPMind1,'uniform',0))==1);
     ind2 = find(cell2mat(cellfun(@(x)sum(ismember(x,i)>0),BPMind2,'uniform',0))==1);
     IND1 = union(ind1,ind2);

     ind1 = find(cell2mat(cellfun(@(x)sum(ismember(x,i)>0),BPM.ind1,'uniform',0))==1);
     ind2 = find(cell2mat(cellfun(@(x)sum(ismember(x,i)>0),BPM.ind2,'uniform',0))==1);
     IND2 = union(ind1,ind2);
     nBPM(i) = length(IND1)/length(IND2);
end
nBPM(find(isnan(nBPM))) = -1;
close all
scatter(-log10(AD_SNP_protective_pv(1,find(nBPM~=-1)))',nBPM(find(nBPM~=-1))');
xlabel('Single SNP association -log10(p-value)')
ylabel('fraction(sig. BPM pv<=0.05 /all BPM)')
saveas(gcf,'scatter_single_SNP_association_BPM_discovery_mhygeSSI_SNPperm.png')

% single SNP association vs BPM
tmp1 = cell2mat(cellfun(@(x)mean(-log10(AD_SNP_protective_pv(1,x))),BPM.ind1,'uniform',0));
tmp2 = cell2mat(cellfun(@(x)mean(-log10(AD_SNP_protective_pv(1,x))),BPM.ind2,'uniform',0));
mean_single_association = max(tmp1,tmp2);
ind = find(BPM.ind1size>=5 & BPM.ind2size>=5);
close all
scatter(mean_single_association(ind),-log10(bpm_local_pv{1}(ind)));
xlabel('average single association')
ylabel('BPM signficance (-log10 p-value)')
saveas(gcf,'scatter_average_single_SNP_association_BPM_discovery_mhygeSSI_SNPperm.png')

ind = find(BPM.ind1size>=5 & BPM.ind2size>=5);
tmp1 = cell2mat(cellfun(@(x)max(-log10(AD_SNP_protective_pv(1,x))),BPM.ind1(ind),'uniform',0));
tmp2 = cell2mat(cellfun(@(x)max(-log10(AD_SNP_protective_pv(1,x))),BPM.ind2(ind),'uniform',0));
max_single_association = max(tmp1,tmp2);
close all
scatter(max_single_association,-log10(bpm_local_pv{1}(ind)));
xlabel('max single association')
ylabel('BPM signficance (-log10 p-value)')
saveas(gcf,'scatter_max_single_SNP_association_BPM_discovery_mhygeSSI.png')

