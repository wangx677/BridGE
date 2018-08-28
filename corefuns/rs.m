function p = rs(M,ind)

% right tail ranksum test between M(:,ind) and M(:, not_ind)
% M is a pxq matrix
% ind is the index for sub matrix
% 

[p,q] = size(M);

[M tieadj]= sparsetiedrank(M(:));
M = reshape(M,[],q);

nx = numel(M(:,ind));
ny = p*q-nx;

srank = M(:,ind);
w = sum(srank(:));
wmean = nx*(nx + ny + 1)/2;
tiescor = 2 * tieadj / ((nx+ny) * (nx+ny-1));
wvar  = nx*ny*((nx + ny + 1) - tiescor)/12;
wc = w - wmean;
z = (wc - 0.5)/sqrt(wvar);
p = normcdf(-z);

