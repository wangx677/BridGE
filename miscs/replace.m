function N = replace(N,O,R)

% FUNCTION N = replace(N,O,R)
% 
% REPLACE replaces all O values in N with corresponding values in R
%
% INPUTS:
%   N - a matrix or vector
%   O - a column vector represents original values
%   R - a column vector represents corresponding replacements for original values
%
%  OUTPUTS:
%   N - new matrix or vector with updated values

[tf,locs] = ismember(N, O);
X = locs(tf);
N(tf) = R(X);
