function removerelatedindividual(genomefile,pi_hat)

data=readtable(genomefile,'filetype','text','Format','%q%q%q%q%q%q%q%q%q%f%q%q%q%q');

data1 = data(:,1:2);
data2 = data(:,3:4);

data2.Properties.VariableNames = data1.Properties.VariableNames;
varname = [data1;data2];
varname = unique(varname);

[tmp idx1] = ismember(data1,varname);
[tmp idx2] = ismember(data2,varname);

v = data.PI_HAT; 
ibd = sparse(idx1,idx2,v,size(varname,1),size(varname,1));
ibd = max(ibd,ibd');

S = sum(ibd>=pi_hat);
S_sorted = sort(S,'descend');

ind2remove = [];
i=1;
while S_sorted(1)~=0
     ind2remove = [ind2remove unique(find(S==S_sorted(1)))];
     ibd(ind2remove,:) = 0;
     ibd(:,ind2remove) = 0;
     S = sum(ibd>=pi_hat);
     S_sorted = sort(S,'descend');
     i = i+1;
end

subject2remove = varname(ind2remove,:);
writetable(subject2remove,'related_subject2remove.txt','Delimiter','\t')     
