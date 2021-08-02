function pv = rankpv(M)

% this script is used to convert each column of M matrix  to a p-value vector 
% based on each member's ranking (high to low) in the vector

M_sorted = sort(M,'descend');
for i=1:size(M,2)
     [~, pv(:,i)] = ismember(M(:,i),M_sorted(:,i),'legacy');
end
pv = pv/size(M,1);
