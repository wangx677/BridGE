function ssM11 = ssm11(data,pheno,alpha1,alpha2,HygeTable)

x = find(pheno==0); % controls
y = find(pheno==1); % cases

data1 = data;
data0 = 1- data;

[p q] = size(data);

clear data

for tt=1:2
     tmp = nchoosek(q,2);
     ssM11{tt} = zeros(1,tmp);
     clear tmp

     a = data1(x,:)'*data1(x,:);
     b = data1(y,:)'*data1(y,:);
     a = triu(a,+1)+tril(a,-1);
     b = triu(b,+1)+tril(b,-1);
     a = squareform(a);
     b = squareform(b);

     if (tt==1)
          ssM11{tt} = computessm(ssM11{tt},1:length(a),x,a,b,HygeTable);
     else
          ssM11{tt} = computessm(ssM11{tt},1:length(a),y,b,a,HygeTable);
     end
     ssM11{tt} = sparse(ssM11{tt});
     ssM11{tt}(isinf(ssM11{tt})) = 16;
end
end

function ssM = computessm(ssM,IND,x,a,b,HygeTable)
% compute snp-snp interactions scores based on HygeTable
test = [ones(length(IND),1)*length(x) a(IND)'+b(IND)' a(IND)'];
[tmp ind]=ismember(test,HygeTable(:,1:3),'rows');
if (nnz(ind==0)>0)
      ind1 = find(ind~=0 & test(:,1)~=0);
      ind = ind(ind1);
      indnew = IND(ind1);
      ssM(indnew) = HygeTable(ind,4);
else
      ssM(IND) = HygeTable(ind,4);
end
end

