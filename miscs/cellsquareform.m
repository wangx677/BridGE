function Z = cellsquareform(Y)

%CELLSQUAREFORM convert a symmetric cell matrix to its vector form
%
% This function has the same purpose as squareform,
% but it can be applied on cell matrix.
% INPUT:
%   Y has to be a symmetric cell matrix, or a diagonal matrix

n = size(Y,1);
Z = Y(tril(true(n),-1));
Z = Z';
