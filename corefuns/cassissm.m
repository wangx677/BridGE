function cassissm(cassifile)

load SNPdataAR.mat

C = strsplit(cassifile,'.');
p = length(SNPdata.rsid);
q = p;

if ismember(C(end),'je')==1
	% joint effect (cassi)
	system(sprintf('grep -v -w NA %s > tmp.je',cassifile)) 
	data = readtable('tmp.je','filetype','text');
     clear tmp*
	system('rm tmp.je')
elseif ismember(C(end),'lr')==1
	% logistic regression (cassi)
	data = readtable(cassifile,'filetype','text');
end

x = data.SNP1;
y = data.SNP2;
z = -log10(data.LR_P);
cassi = sparse(x,y,z,p,q);

log_or = data.LR_LOG_OR;
log_or = sparse(x,y,log_or,p,q);


ssM{1} = cassi.*(log_or<0);
ssM{2} = cassi.*(log_or>0);

clear cassi log_or

for tt=1:2
	ssM{tt} = max(triu(ssM{tt},1),triu(ssM{tt},1)');
     ssM{tt} = squareform(ssM{tt});
     ssM{tt}(isnan(ssM{tt})) = 0;
end


save(sprintf('ssM_%s_%s.%s.mat',C{3},C{1},C{2}),'ssM','-v7.3')
