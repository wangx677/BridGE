function hygessi2pv(ssmFile,n)

% ssmFile: SNP-SNP interaction file name without "_Rxxx.mat"
% n number of permutation

filestr = strsplit(ssmFile,'_');

if strcmp(filestr{end},'combined')~=1
     for R = 0:n
          load(sprintf('%s_R%s.mat',ssmFile,num2str(R)));
%          if exist('ssM_orig','var')~=1
               for tt = 1:2
                    ssMnew{tt}(R+1,:) = ssM{tt};
               end
%          else
%               error('check ssmFile!')
%          end
     end

     for tt=1:2
          ssMpv{tt} = sparse(n+1,size(ssMnew{tt},2));
          for i=1:(n+1)
               ind = find(ssMnew{tt}(i,:)~=0);
               tmp1 = repmat(ssMnew{tt}(i,ind),n+1,1);
               tmp2 = ssMnew{tt}(:,ind) - tmp1;
               ssMpv{tt}(i,ind) = -log10(sum(tmp2>=0)/(n));
               ssMpv{tt}(ssMpv{tt}<0 | ssMnew{tt}<0.1) = 0;
          end          
     end

     clear ssMnew

     for R = 0:n
          load(sprintf('%s_R%s.mat',ssmFile,num2str(R)));
          ssM1 = ssM;
          for tt=1:2
               ssM{tt} = ssMpv{tt}(R+1,:);
          end
          ssM_orig = ssM1;
          clear ssM1
          save(sprintf('%s_R%s.mat',ssmFile,num2str(R)),'ssM','ssM_orig','-v7.3');
     end
else
     for R = 0:n
          file = strjoin(filestr(1:end-1),'_');

          load(sprintf('%s_RR_R%s.mat',file,num2str(R)))
          % check if the ssM interaction score is permutation p-value based
          if exist('ssM_orig','var')==1
               ssM1 = ssM;
          end
          clear ssM ssM_orig

          load(sprintf('%s_DD_R%s.mat',file,num2str(R)))
          if exist('ssM_orig','var')==1
               ssM2 = ssM;
          end
          clear ssM ssM_orig

          load(sprintf('%s_RD_R%s.mat',file,num2str(R)))
          if exist('ssM_orig','var')==1
               ssM3 = ssM;
          end
          clear ssM ssM_orig
     
          if exist('ssM1','var') & exist('ssM2','var') & exist('ssM3','var')
               for tt=1:2
                    maxidx{tt} = zeros(size(ssM1{tt}));
                    ssMtmp{tt} = max(ssM1{tt},ssM2{tt});
                    maxidx{tt}(ssM1{tt} >= ssM2{tt}) = 1;
                    maxidx{tt}(ssM2{tt} > ssM1{tt}) = 2;
                    ssM{tt} = max(ssMtmp{tt},ssM3{tt});
                    maxidx{tt}(ssM3{tt} > ssMtmp{tt}) = 3;
                    maxidx{tt} = maxidx{tt}.*(ssM{tt}>0);
               end
               save(sprintf('%s_combined_R%s.mat',file,num2str(R)),'ssM','maxidx','-append','-v7.3')
               clear ssM maxidx ssMtmp ssM1 ssM2 ssM3
          else
               error('check ssmFile')
          end
     end
end
