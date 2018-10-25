function yeast_flip_positive_negative_info_xls(file)

% this is a temporary script used to flip the 'freq_case' and 'freq_control' labels 
% for files in BPM_WPM_info folder because the case/control setting needs to be flipped

data = readtable(file,'Sheet',1);
data.Properties.VariableNames = {'snps1','genes1','snps2','genes2','GItype','freq_control','freq_case','GI'};
writetable(data,file,'Sheet',1);

tmp = strsplit(file,'_');
if strcmp(tmp{3},'protective')
     tmp{3} = 'risk';
elseif strcmp(tmp{3},'risk')
     tmp{3} = 'protective';
end

newfile = strjoin(tmp,'_');
movefile(file,newfile);
