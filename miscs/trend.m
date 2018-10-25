function [p_CA,stats]=trend(x,m,w)
%COCHRAN_ARM Cochran-Armitage chi-square test for trend in proportions
%   [p_CA,stats]=trend(x,m,w)
%   x - counts matrix, each row presents count for two groups
%           group1        group2
%     -------------------------------
%       x1 x2 ... xm  y1 y2 ... ym     
%
%   m - number of element for additive trend    
%   w - m-by-1 matrix of weights (optional, 1:n by deafult)
%   p_CA - p-value for the main hypothesis
%   stats - structure with overall and Cochran-Armitage chi^2 values, and p
%       values for the test itself (same as p_CA) and for the deviation
%       forom linear trend
%
% The formulas were taken from the books mentioned below and tested on the
% examples from these books, also shown below.
%
% Agresti, Categorical Data Analysis, 2nd ed, pp. 179-182
% x=[ 48  17066;  38  14464;  5  788;  1  126;  1  37 ];
% w=[0 0.5 1.5 4 7];
% 
% Armitage, Berry, Matthews, Statistical Methods in Medical Research, 
% 4th ed, pp. 504-506 
% x=[ 59 97; 10 31; 12 36; 5 28]; w=[0 1 2 3];
% 
% Implemented by Alexander modified by Wen Wang for vectorization
% https://www.mathworks.com/matlabcentral/fileexchange/45186-cochran-armitage-test
% 
 
%% deal with inputs
% check sizes ond orientation of x and w
if nargin<1, error('cochran_arm: Not enough arguments'), end

if size(x,2)~=2*m, x=x'; end;
if size(x,2)~=2*m || size(x,2)<4
    error('cochran_arm: Data must be m x 2')
end

% default: linear weights
if nargin<3, w=1:m; end

if ~isvector(w)
    error('cochran_arm: Weights matrix must be a vector')
end
if size(w,1)~=1, w=w'; end;
if length(w)~=size(x,2)/2
    error('cochran_arm: Sizes of weight matrix and data matrix are incompatible')
end

w = repmat(w,size(x,1),1);

%% computation
% n=sum(x,2); N=sum(n);
% p=x(:,1)./n;
% pp=sum(x(:,1))/N;
% R=sum(x(:,1));
n = [];
for i=1:m
     n = [n sum(x(:,[i,i+m]),2)];
end

N = sum(n,2);
p = x(:,1:m)./n;
pp = sum(x(:,1:m),2)./N;
R = sum(x(:,1:m),2);

%overall chi^2
% x2=1/(pp*(1-pp)) * sum(n.*(p-pp).^2);
x2=1./(pp.*(1-pp)).*sum(n.*(p-repmat(pp,1,m)).^2,2);


%Cochran_armitage chi^2 (aka z^2)
% x2_1_numer=N*(N*sum(x(:,1).*w)-R*sum(n.*w))^2;
% x2_1_denom=R*(N-R)*(N*sum(n.*w.*w)-(sum(n.*w))^2);
% x2_1=x2_1_numer/x2_1_denom;

x2_1_numer=N.*(N.*sum(x(:,1:m).*w,2)-R.*sum(n.*w,2)).^2;
x2_1_denom=R.*(N-R).*(N.*sum(n.*w.*w,2)-(sum(n.*w,2)).^2);
x2_1=x2_1_numer./x2_1_denom;

p_CA=1-chi2cdf(x2_1,1);
p_fit=1-chi2cdf(x2-x2_1,m-2);

%% output
stats.overall_chi2=x2;
stats.cochran_arm_chi2=x2_1;
stats.cochran_arm_p=p_CA;
stats.deviation_from_linear_p=p_fit;
