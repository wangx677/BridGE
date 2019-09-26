function computessi_quant(model,marginal,alpha1,alpha2,plinkCluster2,nWorker,R)

if strcmp(model,'RR')
     load('SNPdataAR.mat')
     for tt=1:2
          ssM{tt} = zeros(size(SNPdata.data,2),size(SNPdata.data,2));
     end     
     tic;
     for i=1:size(SNPdata.data,2)
          for j=i+1:size(SNPdata.data,2)-1
               ind11 = find(SNPdata.data(:,i)==1 & SNPdata.data(:,j) ==1);
               ind10 = find(SNPdata.data(:,i)==1 & SNPdata.data(:,j) ==0);
               ind01 = find(SNPdata.data(:,i)==0 & SNPdata.data(:,j) ==1);
               ind00 = find(SNPdata.data(:,i)==0 & SNPdata.data(:,j) ==0);
               if length(ind11)>=10 & mean(SNPdata.pheno(ind11)) > mean(SNPdata.pheno(ind00)) & mean(SNPdata.pheno(ind11)) > mean(SNPdata.pheno(ind01)) & mean(SNPdata.pheno(ind11)) > mean(SNPdata.pheno(ind10))
                    tmp1 = ranksum(SNPdata.pheno(ind11),SNPdata.pheno(ind00),'tail','right');
                    tmp2 = ranksum(SNPdata.pheno(ind11),SNPdata.pheno(ind10),'tail','right');
                    tmp3 = ranksum(SNPdata.pheno(ind11),SNPdata.pheno(ind01),'tail','right');
                    ssM{1}(i,j) = max(tmp1,max(tmp2,tmp3));
                    ssM{2}(i,j) = max(1-tmp1,max(1-tmp2,1-tmp3));
                end
          end
     end
     toc
end
