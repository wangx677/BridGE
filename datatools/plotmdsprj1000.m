function plotmdsprj1000(mdsFile,popFile,prj1000Pop,outputFile,legendPos)

% Inputs:
% mdsFile: MDS file from plink
% popFile: 1000 project population file (FID, IID, Population)
% prj1000Pop: population that has to be included in the 1000  project 
%             so that we can use to compare with the study population
% outputFile: output file name
% legendPos: legend position
%

% set default legendPos 
if exist('legendPos','var')~=1
     legendPos = 'southeast';
end

% load data
mds = readtable(mdsFile,'filetype','text','Format','%s%s%s%f%f');
popnids = readtable(popFile,'filetype','text','ReadVariableName',0);

% set colors
colorlist = {'m','y','c','r','g','b'};

% map population for project 1000 samples
unique_pop = unique(popnids.Var3);
pop = repmat({'Study'},size(mds,1),1);


j=1;
for i=1:length(unique_pop)
     ind = find(ismember(popnids.Var3,unique_pop{i}));
     ind1 = find(ismember(mds.FID,popnids.Var1(ind)));
     if isempty(ind1)~=1;
           pop(ind1) = unique_pop(i);
           upop{j} = unique_pop{i};
           j = j+1;
     end
end

figure
hold on

ind1 = find(ismember(pop,'Study'));
scatter(mds.C1(ind1),mds.C2(ind1),'+',colorlist{1});

for i=1:length(upop)
     ind1 = find(ismember(pop,upop{i}));
     scatter(mds.C1(ind1),mds.C2(ind1),'+',colorlist{i+1});
end

xlabel('C1');
ylabel('C2');
title('MDS plot')
legend(['Study' upop],'Location',legendPos);
setfig
grid on
grid minor
saveas(gcf,sprintf('%s.pdf',outputFile))
saveas(gcf,sprintf('%s.png',outputFile))


