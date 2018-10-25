function plot_mds(mdsFile,outputFile,groupFile,colorlist)

[FID, IID, SOL, C1, C2] = textread(mdsFile,'%s%s%d%f%f','headerlines',1);

if exist(groupFile,'file')
     [FID1 groups] = textread(groupFile,'%s%s','delimiter','\t');

     for i=1:length(FID)
          ind = find(ismember(FID1,FID{i})==1);
          groupsnew{i} = groups{ind};
     end

     groupsnew(find(ismember(groupsnew,''))) = {'unknown'};
     groups = groupsnew;
end

if exist('groups','var')
     if exist('colorlist','var')~=1
          colorlist={'black','red','green','blue','yellow','#ff1493','#ae017e','#dd3497','#f768a1','#fa9fb5'};
     end
     unique_groups = unique(groups);
     for i=1:length(unique_groups)
          n = nnz(ismember(groups,unique_groups{i}));
          unique_groups_new{i} = sprintf('%s(%d)',unique_groups{i},n);
     end

     if length(unique_groups)<=length(colorlist)
          f = figure()
          hold on
          for i=1:length(unique_groups)
               ind = find(ismember(groups,unique_groups{i}));
               plot(C1(ind),C2(ind),'.','color',colorlist(i))
          end

          setfig
          xlabel('C1')
          ylabel('C2')
          legend(unique_groups_new,'Location','EastOutside','Interpreter', 'none')

     else
          colorlist = dcolors(length(unique_groups));
          f = figure();
          hold on
          for i=1:length(unique_groups)
               ind = find(ismember(groups,unique_groups{i}));
               plot(C1(ind),C2(ind),'.','color',colorlist(i,:))
          end
          
          setfig
          xlabel('C1')
          ylabel('C2')
          legend(unique_groups_new,'Location','EastOutside','Interpreter', 'none')     

     end
else
     plot(C1,C2,'.','color','b');
     setfig
     xlabel('C1')
     ylabel('C2')
end

saveas(gcf,sprintf('%s.pdf',outputFile))
saveas(gcf,sprintf('%s.png',outputFile))
