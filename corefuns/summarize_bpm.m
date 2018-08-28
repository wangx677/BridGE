function summary= summarize_bpm(ssmFile,resultFile,diseasemodel,snppathwayfile)

%% construct subjectxBPM matrix (discovery data and validation data(real and 100 random))
%% 1) subject x snp_pair (show interaction in the BPM) matrix : binary
%% 2) summarize 1) by sum

load(resultFile)
load BPMind.mat
load(ssmFile)
load(snppathwayfile);

load SNPdataAR.mat
trans_data_AR = SNPdata.data;

load SNPdataAD.mat
trans_data_AD = SNPdata.data;

pheno = SNPdata.pheno;

ind = indBPM(find(fdrBPMnew<0.255));

path1 = snpset.pathwaynames([BPM.path1idx BPM.path1idx]);
path2 = snpset.pathwaynames([BPM.path2idx BPM.path2idx]);;

subject_BPM = [];

for i=1:2
	ssM{i} = squareform(ssM{i});
	if strcmp(diseasemodel,'combined')
		maxidx{i} = squareform(maxidx{i});
	elseif strcmp(diseasemodel,'MM')
		maxidx{i} = 2*ones(size(ssM{i}));
	elseif strcmp(diseasemodel,'mm')
		maxidx{i} = ones(size(ssM{i}));
	end
end

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

	   ind1 = BPM.ind1{indtmp};
        ind2 = BPM.ind2{indtmp};
        [m n] = find(ssM{tt}(ind1,ind2)>0.2);

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
	subject_BPM = [subject_BPM output];
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
path1 = path1(ind);
path2 = path2(ind);
FDR = fdrBPMnew';
cases = cases';
controls = controls';
oddsratio = oddsratio';
RorP = RorP';
summary = table(path1, path2, FDR, RorP, cases, controls, oddsratio);

save('summarize_BPM_FDR025.mat','subject_BPM','summary')
