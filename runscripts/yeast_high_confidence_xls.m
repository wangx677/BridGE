function yeast_high_confidence_xls(file)

% this is a temporary script used to flip the 'protective' and 'risk' labels
% for high confidence summary files in /project/csbio/wwang/BridGE/results_collection/yeast
% all because the case/control setting needs to be flipped

data = readtable(file,'Sheet',1);

if isempty(data)~=1;
     data.effect(find(ismember(data.effect,'protective')==1)) = {'tmp'};
     data.effect(find(ismember(data.effect,'risk')==1)) = {'protective'};
     data.effect(find(ismember(data.effect,'tmp')==1)) = {'risk'};
     writetable(data,file,'Sheet',1)
end
