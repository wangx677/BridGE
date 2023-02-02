function cassissm(plinkbimfile,cassifile,gitype,outputfile)

% this function converts CASSI output file to SNP-SNP interaction file that can be used as an input to BridGE
% 
% inputs:
%    cassifile: CASSI output file
%    gitype: genetic interaction type
%            'lin': logistic regression test
%            'lr': logistic regression test 
%            'je': joint effects
%            'awu': adjusted Wu method
%            'afe': adjusted fast epistasis
%            'wz': Wellek Ziegler method
%    outputfile: output file name
%        

% The order of the SNPs is based on the SNP order in the plink file
% Need to know in total how many SNPs are in plink file.
% If there is p SNPs, then the SNP-SNP interaction matrix will be a p x p matrix

% get number of SNPs from plink bim file 
fid=fopen(plinkbimfile);
g = textscan(fid,'%s','delimiter','\n');
fclose(fid)
p = length(g{1});

data = readtable(cassifile,'filetype','text');
x = data.SNP1;
y = data.SNP2;

if strcmp(gitype,'lin')==1 
     % Linear Regression
     z = -log10(data.LIN_P);
     cassi = sparse(x,y,z,p,p);
     lin_beta = data.LIN_BETA;
     lin_beta = sparse(x,y,lin_beta,p,p);
     ssM{1} = cassi.*(lin_beta<0); % protective
     ssM{2} = cassi.*(lin_beta>=0); % risk
elseif strcmp(gitype,'lr')==1
     % logistic regression (cassi)
     z = -log10(data.LR_P);
     cassi = sparse(x,y,z,p,p);
     log_or = data.LR_LOG_OR;
     log_or = sparse(x,y,log_or,p,p);
     ssM{1} = cassi.*(log_or<0); % protective
     ssM{2} = cassi.*(log_or>0); % risk
elseif strcmp(gitype,'je')==1
	% joint effect (cassi)
     z = -log10(data.JE_CC_P);
     cassi = sparse(x,y,z,p,p);
     case_log_or = data.JE_CASE_LOG_OR;
     case_log_or = sparse(x,y,case_log_or,p,p);
     ctrl_log_or = data.JE_CTRL_LOG_OR;
     ctrl_log_or = sparse(x,y,ctrl_log_or,p,p);
     ssM{1} = cassi.*(ctrl_log_or>case_log_or); % protective
     ssM{2} = cassi.*(ctrl_log_or<=case_log_or); % risk
elseif strcmp(gitype,'awu')==1
     % Adjusted Wu
     z = -log10(data.AWU_CC_P);
     cassi = sparse(x,y,z,p,p);
     case_log_or = data.AWU_CASE_LOG_OR;
     case_log_or = sparse(x,y,case_log_or,p,p);
     ctrl_log_or = data.AWU_CTRL_LOG_OR;
     ctrl_log_or = sparse(x,y,ctrl_log_or,p,p);
     ssM{1} = cassi.*(ctrl_log_or>case_log_or); % protective
     ssM{2} = cassi.*(ctrl_log_or<=case_log_or); % risk
elseif strcmp(gitype,'afe')==1
     % Adjusted Fast Epistasis
     z = -log10(data.AFE_CC_P);
     cassi = sparse(x,y,z,p,p);
     case_log_or = data.AFE_CASE_LOG_OR;
     case_log_or = sparse(x,y,case_log_or,p,p);
     ctrl_log_or = data.AFE_CTRL_LOG_OR;
     ctrl_log_or = sparse(x,y,ctrl_log_or,p,p);
     ssM{1} = cassi.*(ctrl_log_or>case_log_or); % protective
     ssM{2} = cassi.*(ctrl_log_or<=case_log_or); % risk
elseif strcmp(gitype,'wz')==1
     % Wellek Ziegler
     z = -log10(data.WZ_CC_P);
     cassi = sparse(x,y,z,p,p);
     case_r = data.WZ_CASE_R;
     case_r = sparse(x,y,case_r,p,p);
     ctrl_r = data.WZ_CTRL_R;
     ctrl_r = sparse(x,y,ctrl_r,p,p);
     ssM{1} = cassi.*(ctrl_r>case_r); % protective
     ssM{2} = cassi.*(ctrl_r<=case_r); % risk
end

for tt=1:2
     ssM{tt}(isnan(ssM{tt})) = 0;
	ssM{tt} = max(ssM{tt},ssM{tt}');
     ssM{tt} = squareform(ssM{tt});
end


save(outputfile,'ssM','-v7.3')
