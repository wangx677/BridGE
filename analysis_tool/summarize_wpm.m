function [subject_WPM subject_WPM_protective subject_WPM_risk summary] = summarize_wpm(ssmFile,BPMindfile,resultFile,diseasemodel,snppathwayfile,fdrcutoff,pcutoff,netcut,outputfile)

%% construct subjectxWPM matrix (discovery data and validation data(real and 100 random))
%% 1) subject x snp_pair (show interaction in the WPM) matrix : binary
%% 2) summarize 1) by sum

load(resultFile)
load(BPMindfile)
load(ssmFile)
load(snppathwayfile);

load SNPdataAR.mat
trans_data_AR = SNPdata.data;

load SNPdataAD.mat
trans_data_AD = SNPdata.data;

pheno = SNPdata.pheno;

ind = find(fdrWPM2<=fdrcutoff & wpm_pv<=pcutoff);
path1 = reshape(WPM.pathway,1,length(WPM.pathway));
path1 = [path1 path1];
path2 = path1;

subject_WPM = [];

for i=1:2
	ssM{i} = squareform(ssM{i});
	if strcmp(diseasemodel,'combined')
		maxidx{i} = squareform(maxidx{i});
	elseif strcmp(diseasemodel,'DD')
		maxidx{i} = 2*ones(size(ssM{i}));
	elseif strcmp(diseasemodel,'RR')
		maxidx{i} = ones(size(ssM{i}));
     elseif strcmp(diseasemodel,'RD')
          maxidx{i} = 3*ones(size(ssM{i}));
	end
end

if isempty(ind)==0
for k = 1:length(ind)
	if ind(k)<=length(path1)/2
		tt = 1;
		indtmp = ind(k);
		RorP{k} = 'protective';
	else
		tt = 2;
		indtmp = ind(k)-length(path1)/2;
		RorP{k} = 'risk';
	end

	   ind1 = WPM.ind{indtmp};
        ind2 = WPM.ind{indtmp};
        
        % for WPM, ssM{tt}(ind1,ind2) is symmetric  
        tmp = tril(ssM{tt}(ind1,ind2),-1);
        [m n] = find(tmp>=netcut);

        ind1 = ind1(m);
        ind2 = ind2(n);

        output = zeros(size(trans_data_AR,1),length(ind1));

        for ii=1:length(ind1)
        	if (maxidx{tt}(ind1(ii),ind2(ii)) == 1)
                	output(:,ii) = trans_data_AR(:,ind1(ii)).* trans_data_AR(:,ind2(ii)) == 1;
                	genotype{k}{ii}='11';
        	elseif (maxidx{tt}(ind1(ii),ind2(ii)) == 2)
                	output(:,ii) = trans_data_AD(:,ind1(ii)).* trans_data_AD(:,ind2(ii)) == 1;
                	genotype{k}{ii}='22';
        	elseif (maxidx{tt}(ind1(ii),ind2(ii)) == 3)
        		tmp1 = trans_data_AD(:,ind1(ii)).* trans_data_AR(:,ind2(ii)) == 1;
                	tmp2 = trans_data_AR(:,ind1(ii)).* trans_data_AD(:,ind2(ii)) == 1;
			if (tt==1)
                		a = nnz(tmp1(pheno==0));
                		hygetmp1 = hygetest(length(pheno),nnz(pheno==0),a,nnz(tmp1));
                		a = nnz(tmp2(pheno==0));
                		hygetmp2 = hygetest(length(pheno),nnz(pheno==0),a,nnz(tmp2));
			else
				a = nnz(tmp1(pheno==1));
                                hygetmp1 = hygetest(length(pheno),nnz(pheno==0),a,nnz(tmp1));
                                a = nnz(tmp2(pheno==1));
                                hygetmp2 = hygetest(length(pheno),nnz(pheno==0),a,nnz(tmp2));
			end

               		if (hygetmp1>hygetmp2)
                        	output(:,ii) = trans_data_AD(:,ind1(ii)).* trans_data_AR(:,ind2(ii)) == 1;
                        	genotype{k}{ii}='21';
                        else
                        	output(:,ii) = trans_data_AR(:,ind1(ii)).* trans_data_AD(:,ind2(ii)) == 1;
                        	genotype{k}{ii}='12';
                	end
   		end
	end
     
	output = sum(output,2);
	subject_WPM = [subject_WPM output];
	cases(k) = nnz(output(pheno==1)>=mean(output));
	controls(k) = nnz(output(pheno==0)>=mean(output));
	if tt==1
		a = controls(k);
		b = cases(k);
		c = nnz(pheno==0)-a;
	        d = nnz(pheno==1)-b;
	elseif tt==2
		a = cases(k);
                b = controls(k);
		c = nnz(pheno==1)-a;
	        d = nnz(pheno==0)-b;
        end

	oddsratio(k) = (a/c)/(b/d);

end
path = reshape(path1(ind),length(ind),1);
FDR = reshape(fdrWPM2(ind),length(ind),1);
permP = reshape(wpm_pv(ind),length(ind),1);
cases = reshape(cases,length(ind),1)/nnz(pheno==1);
controls = reshape(controls,length(ind),1)/nnz(pheno==0);
oddsratio = reshape(oddsratio,length(ind),1);
RorP = reshape(RorP,length(ind),1);
summary = table(path, FDR, permP, RorP, cases, controls, oddsratio);
[tmp1 tmp2] = sort(permP,'ascend');
summary = summary(tmp2,:);


pheno = SNPdata.pheno;
pheno = table(pheno);
subject_WPM = array2table(subject_WPM);
subject_WPM.Properties.VariableNames = arrayfun(@(x)sprintf('WPM%s',num2str(x)),1:length(ind),'uniform',0);

subject_WPM_protective = subject_WPM(:,find(ismember(RorP,'protective')));
subject_WPM_risk = subject_WPM(:,find(ismember(RorP,'risk')));

subject_WPM = [pheno subject_WPM];
subject_WPM_protective = [pheno subject_WPM_protective];
subject_WPM_risk = [pheno subject_WPM_risk];

if exist('outputfile','var')
     save(outputfile,'subject_WPM','summary')
end

else
     subject_WPM_protective = [];
     subject_WPM_risk = [];
     subject_WPM = [];
     summary = [];
end
