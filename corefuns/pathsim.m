function sim = pathsim(pathind)

% FUNCTION sim = pathsim(pathind)
%
% PATHSIM computes pairwise similarity (overlap coefficient)
% for given set of indexes

sim = zeros(length(pathind),length(pathind));

for i=1:length(pathind)
     for j=i+1:length(pathind)
	     sim(i,j) = length(intersect(pathind{i},pathind{j}))/min(length(pathind{i}),length(pathind{j}));
	end
end

sim = max(sim,sim');
end
