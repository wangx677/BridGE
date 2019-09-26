function M = repmatcell(M,X,Y)

[m n] = size(M);

mm = M;

for x=1:X-1
	M = [M;mm];
end

mm = M;
for y=1:Y-1
	M = [M mm];
end
	

