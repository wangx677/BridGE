function yeast_flip_positive_negative_xls(file)

% this is a temporary script used to flip the 'protective' and 'risk' labels
% for files in all output*.xls files because the case/control setting needs to be flipped

try
     data = readtable(file,'Sheet',3);
catch
     data = [];
end

if isempty(data)~=1;
     data.eff_bpm(find(ismember(data.eff_bpm,'protective')==1)) = {'tmp'};
     data.eff_bpm(find(ismember(data.eff_bpm,'risk')==1)) = {'protective'};
     data.eff_bpm(find(ismember(data.eff_bpm,'tmp')==1)) = {'risk'};
     writetable(data,file,'Sheet',3)
end

try
     data = readtable(file,'Sheet',4);
catch
     data = [];
end

if isempty(data)~=1;
     data.eff_wpm(find(ismember(data.eff_wpm,'protective')==1)) = {'tmp'};
     data.eff_wpm(find(ismember(data.eff_wpm,'risk')==1)) = {'protective'};
     data.eff_wpm(find(ismember(data.eff_wpm,'tmp')==1)) = {'risk'};
     writetable(data,file,'Sheet',4);
end

try
     data = readtable(file,'Sheet',5);
catch
     data = [];
end

if isempty(data)~=1;
     data.eff_path(find(ismember(data.eff_path,'protective')==1)) = {'tmp'};
     data.eff_path(find(ismember(data.eff_path,'risk')==1)) = {'protective'};
     data.eff_path(find(ismember(data.eff_path,'tmp')==1)) = {'risk'};
     writetable(data,file,'Sheet',5)
end

