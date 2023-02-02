projectdirs={'project_ALS_NINDS_phs000101_OmniExpress_GOIset1','project_ALS_NINDS_phs000101_HumanHap_GOIset1','project_ALS_Finland_phs000344_GOIset1','project_ALS_Irish_phs000127_GOIset1'};
fdrcutoff=0.25;
outputdir='/project/csbio/wwang/BridGE/results_collection/ALS';

% validcutoff=1:(length(projectdirs)-1);
% tail={'one','two'};
% GI={'mhygeSSI','hygeSSI'};
% fdrcutoff = 0.05:0.05:0.25;
m=1;

for fdrcutoff = 0.05:0.05:0.25
     for GI = {'mhygeSSI','hygeSSI'};
          for tail = {'one','two'}
               for validcutoff=1:(length(projectdirs)-1);
                    for R=0:1000
                    [number_GI_risk(R+1),  number_GI_validated_risk(R+1), number_GI_protective(R+1), number_GI_validated_protective(R+1)] = summarize_result_gene(projectdirs,GI{1},R,fdrcutoff,outputdir,validcutoff,tail{1});
                         if R==0 & number_GI_validated_risk(1) ==0 & number_GI_validated_protective(1)==0
                              break
                         end
                    end
                    GI_all(m) = GI;
                    tail_all(m) = tail;
                    validcutoff_all(m) = validcutoff;
                    fdrcutoff_all(m) = fdrcutoff;

                    number_GI_discovered_risk_all(m) = number_GI_risk(1);
                    number_GI_validated_risk_all(m) = number_GI_validated_risk(1);
                    if number_GI_discovered_risk_all(m)>0
                         pv_risk(m) = nnz(number_GI_validated_risk(2:end)>=number_GI_validated_risk(1))/(length(number_GI_validated_risk)-1);
                    else
                         pv_risk(m)=1;
                    end

                    number_GI_discovered_protective_all(m) = number_GI_protective(1);
                    number_GI_validated_protective_all(m) = number_GI_validated_protective(1);
                    if number_GI_validated_protective_all(m)>0
                         pv_protective(m) = nnz(number_GI_validated_protective(2:end)>=number_GI_validated_protective(1))/(length(number_GI_validated_protective)-1);              
                    else
                         pv_protective(m) = 1;
                    end

                    number_GI_discovered_combined(m) = number_GI_risk(1) + number_GI_protective(1);
                    number_GI_validated_combined_all(m) = number_GI_validated_risk(1)+number_GI_validated_protective(1);
                    if number_GI_validated_combined_all(m)>0
                         pv_combined(m) = nnz((number_GI_validated_risk(2:end)+number_GI_validated_protective(2:end))>=(number_GI_validated_risk(1)+number_GI_validated_protective(1)))/(length(number_GI_validated_risk)-1);
                    else
                          pv_combined(m)=1;
                    end
                    clear number_GI_risk number_GI_validated_risk number_GI_protective number_GI_validated_protective
                    m = m+1
               end
          end
     end
end


GI = reshape(GI_all,length(GI_all),1);
fdrcutoff = reshape(fdrcutoff_all,length(fdrcutoff_all),1);
tail = reshape(tail_all,length(tail_all),1);
validcutoff = reshape(validcutoff_all,length(validcutoff_all),1);

NO_GI_risk_discover = reshape(number_GI_discovered_risk_all,length(number_GI_discovered_risk_all),1);
NO_GI_risk_valid = reshape(number_GI_validated_risk_all,length(number_GI_validated_risk_all),1);
pv_risk = reshape(pv_risk,length(pv_risk),1);

NO_GI_protective_discover = reshape(number_GI_discovered_protective_all,length(number_GI_discovered_protective_all),1);
NO_GI_protective_valid = reshape(number_GI_validated_protective_all,length(number_GI_validated_protective_all),1);
pv_protective = reshape(pv_protective,length(pv_protective),1);

NO_GI_combined_disocver = reshape(number_GI_discovered_combined,length(number_GI_discovered_combined),1);
NO_GI_combined_valid = reshape(number_GI_validated_combined_all,length(number_GI_validated_combined_all),1);
pv_combined = reshape(pv_combined,length(pv_combined),1);


output = table(GI,fdrcutoff,validcutoff,tail,NO_GI_protective_discover,NO_GI_protective_valid,pv_protective,NO_GI_risk_discover,NO_GI_risk_valid,pv_risk,NO_GI_combined_disocver,NO_GI_combined_valid,pv_combined);
writetable(output,'/project/csbio/wwang/BridGE/results_collection/ALS/results_GOI_summary.xls')

