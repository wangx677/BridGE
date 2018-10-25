function [r, tieadj] = sparsetiedrank(x)

% this is a function for sparse data which includes a majority of zeros and positive values
N = length(x);
x = reshape(x,N,1);

% [sx, rowidx]= sort(x);
ind0 = find(x==0);
ind1 = find(x~=0);
sx0 = zeros(length(ind0),1);
[sx1 rowidx1] = sort(x(ind1));
sx = [sx0; sx1];
rowidx = [ind0; ind1(rowidx1)];

ranks = 1:N;

ties = sx(1:N-1) >= sx(2:N);
tieloc = find(ties);

nottie = find(~ties);
nottie = nottie(find(nottie>tieloc(end)));
if isempty(nottie)
     nottie = N;
else
     nottie = nottie(1);
end

locstart = [tieloc(1); tieloc(find(tieloc(2:end)-tieloc(1:end-1)>1)'+1)];
locend = [tieloc(find(tieloc(2:end)-tieloc(1:end-1)>1)')+1;nottie];

ntied = locend-locstart+1;
tieadj = sum(ntied.*(ntied-1).*(ntied+1)/2);
ranktmp = ((locstart+locend).*(locend-locstart+1)/2)./(locend-locstart+1);

for i=1:length(locstart)
     ranks(locstart(i):locend(i)) = ranktmp(i);
end
r(rowidx) = ranks;
r = reshape(r,size(x));
