model = {'RR','DD','RD','combined'};

for j=1:length(model)
     for i=0:1000
          load(sprintf('ssM_mhygeSSI_alpha10.05_alpha20.05_%s_R%s.mat',model{j},num2str(i)));
          for tt=1:2
               if j==1
                    pbody_density_RR{tt}(i+1) = mean(ssM{tt});
               elseif j==2
                    pbody_density_DD{tt}(i+1) = mean(ssM{tt});
               elseif j==3
                    pbody_density_RD{tt}(i+1) = mean(ssM{tt});
               elseif j==4
                    pbody_density_combined{tt}(i+1) = mean(ssM{tt});
               end
          end
     end
end

for tt=1:2
     pbody_density_RR_pv{tt} = nnz(pbody_density_RR{tt}(2:end)>=pbody_density_RR{tt}(1))/1000;
     pbody_density_DD_pv{tt} = nnz(pbody_density_DD{tt}(2:end)>=pbody_density_DD{tt}(1))/1000; 
     pbody_density_RD_pv{tt} = nnz(pbody_density_RD{tt}(2:end)>=pbody_density_RD{tt}(1))/1000;
     pbody_density_combined_pv{tt} = nnz(pbody_density_combined{tt}(2:end)>=pbody_density_combined{tt}(1))/1000;
end

save('results_pbody_WPM.mat','pbody_density_*_pv')
