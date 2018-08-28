function [ measures ] = chi2test(tables)

% FUNCTION [ measures ] = chi2(tables)
% CHI2TEST  compute chi-square from 2 by 2 contigency tables 
%        
% The contingency table looks as follows:
%
%                   q (cases)       not q (controls)  row totals
%                 ----------------------------------------------
%            p     f11 (a)          f10 (b)            f1plus
%        not p     f01 (c)          f00 (d)            f0plus
%                 ----------------------------------------------
% column totals   fplus1 (cases)    fplus0 (controls)    N
%
% In terms of other common notation, a = f11, b = f10, c = f01, and d = f00.
% It is assumed that tables is a four by number of tables arrays, is passed 
% in the order f11, f10, f01, f00 or a b c d.
% 
num_tables = size( tables, 1 );
measures = zeros( num_tables, 1);

f11 = tables( :, 1 );
f10 = tables( :, 2 );
f01 = tables( :, 3 );
f00 = tables( :, 4 );

f1plus = f11 + f10;
f0plus = f01 + f00;
fplus1 = f11 + f01;
fplus0 = f10 + f00;
N = f1plus + f0plus;


thesig = 0;
xx = sum(tables(1,:));
if xx < 50
      return;
end

c1 = fplus1.*f1plus/xx;
c2 = fplus1.*f0plus/xx;
c3 = fplus0.*f1plus/xx;
c4 = fplus0.*f0plus/xx;

five = 5;
c = (c1>=five).*(c2>=five).*(c3>=five).*(c4>=five);

u1 = ((f11).* (f00) - (f10) .* (f01)).^2;
u2 = (f11) + (f00) + (f10) + (f01);
d1 = (f11) + (f10);
d2 = (f11) + (f01);
d3 = (f00) + (f10);
d4 = (f00) + (f01);
measures = (u1.*u2)./(d1.*d2.*d3.*d4);

measures(isnan(measures))=0;

