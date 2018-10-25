load project_MSSNG_hiseqx_CP/results_ssM_hygeSSI_alpha10.05_alpha20.05_combined_R0.mat

fdr_bpm_combined = fdrBPM2;

fdr_wpm_combined = fdrWPM2;

fdr_path_combined = fdrPATH2;

clearvars -except fdr_*

load project_MSSNG_hiseqx_CP/BPMind.mat
load project_MSSNG_hiseqx_CP/snp_pathway_min10_max300.mat

ind_bpm = find(fdr_bpm_combined<=0.25);
ind_wpm = find(fdr_wpm_combined<=0.3);
ind_path = find(fdr_path_combined<=0.25)

path1 = [WPM.pathway(BPM.path1idx) WPM.pathway(BPM.path1idx)];
path2 = [WPM.pathway(BPM.path2idx) WPM.pathway(BPM.path2idx)];

bpms = [path1(ind_bpm)' path2(ind_bpm)'];

path = [WPM.pathway WPM.pathway];
wpms = [path(ind_wpm)' path(ind_wpm)']

network = [bpms;wpms];

% check redundancey

overlap1 = zeros(length(ind_bpm),length(ind_bpm));
overlap2 = zeros(length(ind_bpm),length(ind_bpm));
bpmind1 = [BPM.ind1 BPM.ind1];
bpmind2 = [BPM.ind2 BPM.ind2];

for i=1:length(ind_bpm)
     for j=i+1:length(ind_bpm)-1
          a = length(intersect(bpmind1{ind_bpm(i)},bpmind1{ind_bpm(j)}));
          b = length(intersect(bpmind2{ind_bpm(i)},bpmind2{ind_bpm(j)}));
          n1 = a*b;
          a = length(intersect(bpmind1{ind_bpm(i)},bpmind2{ind_bpm(j)}));
          b = length(intersect(bpmind2{ind_bpm(i)},bpmind1{ind_bpm(j)}));
          n2 = a*b;
          n = min(length(bpmind1{ind_bpm(i)})*length(bpmind2{ind_bpm(i)}), length(bpmind1{ind_bpm(j)})*length(bpmind2{ind_bpm(j)}));
          overlap1(i,j) = n1/n;
          overlap2(i,j) = n2/n;
     end
end

overlap_bpm = max(overlap1,overlap2);
overlap_bpm = max(overlap_bpm,overlap_bpm');
clear overlap1 overlap2
overlap_bpm = overlap_bpm >0.5;

[noRD,group_bpm] = graphconncomp(sparse(overlap_bpm));
group_bpm = group_bpm';

wpmind = [WPM.ind WPM.ind];
overlap1 = zeros(length(ind_wpm),length(ind_wpm));

for i=1:length(ind_wpm)
     for j=i+1:length(ind_wpm)-1
          a = length(intersect(wpmind{ind_wpm(i)},wpmind{ind_wpm(j)}));
          n1 = a*a;
          n = length(wpmind{ind_wpm(i)})*length(wpmind{ind_wpm(j)});
          overlap1(i,j) = n1/n;
     end
end

overlap_wpm = max(overlap1,overlap1');
clear overlap1 
overlap_wpm = overlap_wpm >0.5;

[noRD,group_wpm] = graphconncomp(sparse(overlap_wpm));
group_wpm = group_wpm'+max(group_bpm);


group = [group_bpm; group_wpm];

pathway1 = network(:,1);
pathway2 = network(:,2);

j=0;
for i=1:length(ind_bpm)
     if ind_bpm(i)>length(BPM.size)
          j = j+1;
          RorP{j} = 'risk';
     else
          j = j+1;
          RorP{j} = 'protective';
     end
end

for i=1:length(ind_wpm)
     if ind_wpm(i) > length(WPM.size)
          j = j+1;
          RorP{j} = 'risk';
     else
          j = j+1;
          RorP{j} = 'protective';
     end
end

RorP = RorP';
     
network = table(pathway1, pathway2, RorP, group);
writetable(network,'MSSNG_Simon_Grant_network.csv')


pathways = unique([bpms; wpms]);
for i=1:length(pathways)
     ind_pathways(i) = find(ismember(WPM.pathway,pathways{i}));
end

wpmind = [WPM.ind];
overlap1 = zeros(length(ind_pathways),length(ind_pathways));

for i=1:length(ind_pathways)
     for j=i+1:length(ind_pathways)-1
          a = length(intersect(wpmind{ind_pathways(i)},wpmind{ind_pathways(j)}));
          n1 = a*a;
          n = length(wpmind{ind_pathways(i)})*length(wpmind{ind_pathways(j)});
          overlap1(i,j) = n1/n;
     end
end

overlap_pathways = max(overlap1,overlap1');
clear overlap1
overlap_pathways = overlap_pathways >0.5;

[noRD,group_pathways] = graphconncomp(sparse(overlap_pathways));

group_pathways = group_pathways';

save collect_MSSNG_Simon_Grant network group_pathways pathways
