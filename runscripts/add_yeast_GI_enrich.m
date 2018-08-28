function add_yeast_GI_enrich(file)

enrich_EE = readtable('/project/csbio/wwang/JHOU_Yeast/Data File S13_Genetic interaction enrichment among protein complexes.xlsx','Sheet','ExE complex network_all');
enrich_EN = readtable('/project/csbio/wwang/JHOU_Yeast/Data File S13_Genetic interaction enrichment among protein complexes.xlsx','Sheet','ExN complex network_all');
enrich_NN = readtable('/project/csbio/wwang/JHOU_Yeast/Data File S13_Genetic interaction enrichment among protein complexes.xlsx','Sheet','NxN complex network_all');

     
     try
          data_bpm = readtable(file,'Sheet','Sheet3');
     catch
          data_bpm = '';
     end

     try 
          data_wpm = readtable(file,'Sheet','Sheet4');
     catch
          data_wpm = '';
     end
         

     if isempty(data_bpm)~=1 
     for i=1:size(data_bpm,1)
          ind1 = find(ismember(enrich_EE.QueryComplex,data_bpm.path1{i}) & ismember(enrich_EE.ArrayComplex,data_bpm.path2{i}));
          ind2 = find(ismember(enrich_EE.QueryComplex,data_bpm.path2{i}) & ismember(enrich_EE.ArrayComplex,data_bpm.path1{i}));
          if (isempty(ind1)~=1 | isempty(ind2)~=1)
          GI_enrich_pv_EE_neg(i) = min([enrich_EE.NegativeInteractionEnrichment_Pvalue_(ind1),enrich_EE.NegativeInteractionEnrichment_Pvalue_(ind2)]);
          GI_enrich_pv_EE_pos(i) = min([enrich_EE.PositiveInteractionEnrichment_Pvalue_(ind1),enrich_EE.PositiveInteractionEnrichment_Pvalue_(ind2)]);
          else
          GI_enrich_pv_EE_neg(i) = nan;
          GI_enrich_pv_EE_pos(i) = nan;
          end

          ind1 = find(ismember(enrich_EN.QueryComplex,data_bpm.path1{i}) & ismember(enrich_EN.ArrayComplex,data_bpm.path2{i}));
         ind2 = find(ismember(enrich_EN.QueryComplex,data_bpm.path2{i}) & ismember(enrich_EN.ArrayComplex,data_bpm.path1{i}));
          if (isempty(ind1)~=1 | isempty(ind2)~=1)
          GI_enrich_pv_EN_neg(i) = min([enrich_EN.NegativeInteractionEnrichment_Pvalue_(ind1),enrich_EN.NegativeInteractionEnrichment_Pvalue_(ind2)]);
          GI_enrich_pv_EN_pos(i) = min([enrich_EN.PositiveInteractionEnrichment_Pvalue_(ind1),enrich_EN.PositiveInteractionEnrichment_Pvalue_(ind2)]);
          else
          GI_enrich_pv_EN_neg(i) = nan;
          GI_enrich_pv_EN_pos(i) = nan;
          end

          ind1 = find(ismember(enrich_NN.QueryComplex,data_bpm.path1{i}) & ismember(enrich_NN.ArrayComplex,data_bpm.path2{i}));
          ind2 = find(ismember(enrich_NN.QueryComplex,data_bpm.path2{i}) & ismember(enrich_NN.ArrayComplex,data_bpm.path1{i}));

          if (isempty(ind1)~=1 | isempty(ind2)~=1)
          GI_enrich_pv_NN_neg(i) = min([enrich_NN.NegativeInteractionEnrichment_Pvalue_(ind1),enrich_NN.NegativeInteractionEnrichment_Pvalue_(ind2)]);
          GI_enrich_pv_NN_pos(i) = min([enrich_NN.PositiveInteractionEnrichment_Pvalue_(ind1),enrich_NN.PositiveInteractionEnrichment_Pvalue_(ind2)]);
          else
          GI_enrich_pv_NN_neg(i) = nan;
          GI_enrich_pv_NN_pos(i) = nan;
          end
     end
     
     GI_enrich_pv_EE_neg = reshape(GI_enrich_pv_EE_neg,length(GI_enrich_pv_EE_neg),1);
     GI_enrich_pv_EE_pos = reshape(GI_enrich_pv_EE_pos,length(GI_enrich_pv_EE_pos),1);
     GI_enrich_pv_EN_neg = reshape(GI_enrich_pv_EN_neg,length(GI_enrich_pv_EN_neg),1);
     GI_enrich_pv_EN_pos = reshape(GI_enrich_pv_EN_pos,length(GI_enrich_pv_EN_pos),1);
     GI_enrich_pv_NN_neg = reshape(GI_enrich_pv_NN_neg,length(GI_enrich_pv_NN_neg),1);
     GI_enrich_pv_NN_pos = reshape(GI_enrich_pv_NN_pos,length(GI_enrich_pv_NN_pos),1);

     bpm_enrich = table(GI_enrich_pv_EE_neg,GI_enrich_pv_EE_pos,GI_enrich_pv_EN_neg,GI_enrich_pv_EN_pos,GI_enrich_pv_NN_neg,GI_enrich_pv_NN_pos);
     
     data_bpm_new = [data_bpm(:,1:end-2) bpm_enrich data_bpm(:,end-1:end)];
     writetable(data_bpm_new,file,'Sheet','Sheet3');

     clear ind1 ind2 GI_enrich_pv_EE_neg GI_enrich_pv_EE_pos GI_enrich_pv_EN_neg GI_enrich_pv_EN_pos GI_enrich_pv_NN_neg GI_enrich_pv_NN_pos
     end


     if isempty(data_wpm)~=1
     for i=1:size(data_wpm,1)
          ind = find(ismember(enrich_EE.QueryComplex,data_wpm.path{i}) & ismember(enrich_EE.ArrayComplex,data_wpm.path{i}));
         
          if isempty(ind)~=1
          GI_enrich_pv_EE_neg(i) = min(enrich_EE.NegativeInteractionEnrichment_Pvalue_(ind));
          GI_enrich_pv_EE_pos(i) = min(enrich_EE.PositiveInteractionEnrichment_Pvalue_(ind));
          else
          GI_enrich_pv_EE_neg(i) = nan;
          GI_enrich_pv_EE_pos(i) = nan;
          end

          ind = find(ismember(enrich_EN.QueryComplex,data_wpm.path{i}) & ismember(enrich_EN.ArrayComplex,data_wpm.path{i}));

          if isempty(ind)~=1
          GI_enrich_pv_EN_neg(i) = min(enrich_EN.NegativeInteractionEnrichment_Pvalue_(ind));
          GI_enrich_pv_EN_pos(i) = min(enrich_EN.PositiveInteractionEnrichment_Pvalue_(ind));
          else
          GI_enrich_pv_EN_neg(i) = nan;
          GI_enrich_pv_EN_pos(i) = nan;
          end

          ind = find(ismember(enrich_NN.QueryComplex,data_wpm.path{i}) & ismember(enrich_NN.ArrayComplex,data_wpm.path{i}));

          if isempty(ind)~=1
          GI_enrich_pv_NN_neg(i) = min(enrich_NN.NegativeInteractionEnrichment_Pvalue_(ind));
          GI_enrich_pv_NN_pos(i) = min(enrich_NN.PositiveInteractionEnrichment_Pvalue_(ind));
          else
          GI_enrich_pv_NN_neg(i) = nan;
          GI_enrich_pv_NN_pos(i) = nan;
          end
     end
     
     GI_enrich_pv_EE_neg = reshape(GI_enrich_pv_EE_neg,length(GI_enrich_pv_EE_neg),1);
     GI_enrich_pv_EE_pos = reshape(GI_enrich_pv_EE_pos,length(GI_enrich_pv_EE_pos),1);
     GI_enrich_pv_EN_neg = reshape(GI_enrich_pv_EN_neg,length(GI_enrich_pv_EN_neg),1);
     GI_enrich_pv_EN_pos = reshape(GI_enrich_pv_EN_pos,length(GI_enrich_pv_EN_pos),1);
     GI_enrich_pv_NN_neg = reshape(GI_enrich_pv_NN_neg,length(GI_enrich_pv_NN_neg),1);
     GI_enrich_pv_NN_pos = reshape(GI_enrich_pv_NN_pos,length(GI_enrich_pv_NN_pos),1);

     wpm_enrich = table(GI_enrich_pv_EE_neg,GI_enrich_pv_EE_pos,GI_enrich_pv_EN_neg,GI_enrich_pv_EN_pos,GI_enrich_pv_NN_neg,GI_enrich_pv_NN_pos);

     data_wpm_new = [data_wpm(:,1:end-1) wpm_enrich data_wpm(:,end)];
     writetable(data_wpm_new,file,'Sheet','Sheet4');

     clear ind1 ind2 GI_enrich_pv_EE_neg GI_enrich_pv_EE_pos GI_enrich_pv_EN_neg GI_enrich_pv_EN_pos GI_enrich_pv_NN_neg GI_enrich_pv_NN_pos
     end
