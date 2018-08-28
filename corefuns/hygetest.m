function [logpv pv] = hygetest(n,d,k,m)

% FUNCTION [logpv pv] = hygetest(n,d,k,m)
%
% HYGETEST computes Hypergeometric cumulative distribution for
% the given inputs.
%
% INPUTS:
% n - is the population size
% d - is the number of draws
% k - is the number of successes
% m - is the number of success states in the population
%
% OUTPUTS:
% logpv - negative log10(p-value)
% pv - p-value

pv = 1-hygecdf(k-1,n,m,d);
logpv = -log10(pv);
