function plot_BPM_combined(pathname1,pathname2,model,riskornot,outputname,snppathwayfile,BPMindfile,denscut)

% pathname1='GPI-anchor transamidase complex';
% pathname2='nuclear pore';
% model='DD';
% riskornot=1;
% outputname='GPI_nuclear_pore.pdf'
% snppathwayfile='snp_pathway_min5_max300.mat';
% BPMindfile='BPMind.mat'
% denscut=0.2;

load(sprintf('ssM_hygeSSI_alpha10.05_alpha20.05_%s_R0.mat',model));

if (riskornot==0)
	MM = squareform(ssM{1}>denscut);
     if isequal(model,'combined')
          MMmax = squareform(maxidx{1});
     end
elseif (riskornot==1)
	MM = squareform(ssM{2}>denscut);
     if isequal(model,'combined')
          MMmax = squareform(maxidx{2});
     end
end

load(snppathwayfile)
load(BPMindfile)

load SNPdataAR.mat
data_AR = SNPdata.data;
load SNPdataAD.mat
data_AD = SNPdata.data;
pheno = SNPdata.pheno;
chromoID = SNPdata.chr;

ii = find(ismember(snpset.pathwaynames,pathname1)==1);
jj = find(ismember(snpset.pathwaynames,pathname2)==1);
A = zeros(length(snpset.pathwaynames),length(snpset.pathwaynames));
A(ii,jj) = 1;
A(jj,ii) = 1;
A = squareform(A);
nn = find(A==1);

ind1 = BPM.ind1{nn};
ind2 = BPM.ind2{nn};

clear A nn

d1 = sum2(MM(ind1,:))/numel(MM(ind1,:));
d2 = sum2(MM(ind2,:))/numel(MM(ind2,:));
BPM_expected = max(d1,d2);
BPM_observed = sum2(MM(ind1,ind2))/numel(MM(ind1,ind2));

BPM_enrich_pv = max(rs(MM(ind1,:),ind2),rs(MM(ind2,:),ind1));

if (length(ind1)>length(ind2))
    tmp1 = ind1;
    tmp2 = ind2;
    ind1 = tmp2;
    ind2 = tmp1;
    tmp1 = pathname1;
    tmp2 = pathname2;
    pathname1 = tmp2;
    pathname2 = tmp1;
    clear tmp1 tmp2
end

a=sum(data_AR(pheno==1,:));
b=sum(data_AR(pheno==0,:));
c=nnz(pheno==1)-a;
d=nnz(pheno==0)-b;
hyge=hygetest(length(pheno),nnz(pheno==1),a,a+b);
hyge_tmp1=hyge;

a=sum(data_AD(pheno==1,:));
b=sum(data_AD(pheno==0,:));
c=nnz(pheno==1)-a;
d=nnz(pheno==0)-b;
hyge=hygetest(length(pheno),nnz(pheno==1),a,a+b);
hyge_tmp2=hyge;

MM = MM(ind1, ind2);
if exist('MMmax','var')
     MMmax = MMmax(ind1,ind2);
     MMmax(MM==0) = 0;
end

hyge_ind1_mm = hyge_tmp1(ind1);
hyge_ind1_MM = hyge_tmp2(ind1);

hyge_ind2_mm = hyge_tmp1(ind2);
hyge_ind2_MM = hyge_tmp2(ind2);

chromoID_ind1 = chromoID(ind1);
chromoID_ind2 = chromoID(ind2);

if exist('MMmax','var')
     test_ind1 = arrayfun(@(x)unique(MMmax(x,:)),1:length(ind1),'uniform',0);
     test_ind2 = arrayfun(@(x)unique(MMmax(:,x)),1:length(ind2),'uniform',0);
end

if strcmp(model,'combined')
     for i=1:length(test_ind1)
	     tmp = test_ind1{i}(2:end);
	     if length(tmp)>0
	          for j=1:length(tmp);
	     	     if (tmp(j)==1)
		     	     tmp_p(j) = hyge_ind1_mm(i);
		          elseif (tmp(j)==2)
			          tmp_p(j) = hyge_ind1_MM(i);
		          elseif (tmp(j)==3)
			          tmp_p(j) = max(hyge_ind1_mm(i),hyge_ind1_MM(i));
		          end
	          end
	          hyge_ind1(i) = max(tmp_p(j));
	          clear tmp_p
	     end
     end
elseif strcmp(model,'RR')
     hyge_ind1 = hyge_ind1_mm;
elseif strcmp(model,'DD')
     hyge_ind1 = hyge_ind1_MM;
end
 
if strcmp(model,'combined')
     for i=1:length(test_ind2)
          tmp = test_ind2{i}(2:end);
          if length(tmp)>0
               for j=1:length(tmp);
                    if (tmp(j)==1)
                         tmp_p(j) = hyge_ind2_mm(i);
                    elseif (tmp(j)==2)
                         tmp_p(j) = hyge_ind2_MM(i);
                    elseif (tmp(j)==3)
                         tmp_p(j) = max(hyge_ind2_mm(i),hyge_ind2_MM(i));
                    end
               end
               hyge_ind2(i) = max(tmp_p(j));
               clear tmp_p
           end
     end
elseif strcmp(model,'RR')
     hyge_ind2 = hyge_ind1_mm;
elseif strcmp(model,'DD')
     hyge_ind2 = hyge_ind2_MM;
end


hyge_ind1(hyge_ind1<=0.01) = 0.01;
hyge_ind2(hyge_ind2<=0.01) = 0.01; 

rgb22={'Violet','Turquoise','YellowGreen','Khaki','Salmon','HotPink','Goldenrod',...
    'SteelBlue','LightSeaGreen','Peru','IndianRed','LightSlateGray','DarkOrchid','CadetBlue',...
    'SeaGreen','Moccasin','DimGray','MediumPurple','DodgerBlue','DarkKhaki','RosyBrown','Darkolivegreen'};


for i=22:-1:1
    rgbcc(i,:) = rgb(rgb22{i});
end


f = figure;
hold on
set(0,'DefaultAxesColorOrder',rgbcc)
tmp = rand(22,2);
hdl = plot(tmp','LineWidth',1);
%h = gridLegend(hdl,6,{'Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7',...
%     'Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15',...
%     'Chr16','Chr17','Chr18','Chr19','Chr20','Chr21','Chr22'},'color','none')
%
%set(h,'box','off')
%lh = findobj(h,'type','line');
%set(lh,'linewidth',4);
delete(hdl)

%p1 = 0;p2 = 2;
p1 = 0;p2 = 1;

idx1 = find(sum(MM,2)~=0);
idx2 = find(sum(MM,1)~=0);
MM = MM(idx1,idx2);
idx1 = 1:length(idx1);
idx2 = 1:length(idx2);
hyge_ind1 = hyge_ind1(idx1);
hyge_ind2 = hyge_ind2(idx2);
chromoID_ind1 = chromoID_ind1(idx1);
chromoID_ind2 = chromoID_ind2(idx2);


bb = (max(length(idx1),length(idx2))-min(length(idx1),length(idx2)))/2;
max1 = 0;
for i=1:22
    tt = find(chromoID_ind2(idx2)==i);
    if (nnz(tt)~=0)
    line([p2*ones(size(idx2(tt)));hyge_ind2(idx2(tt))/2+p2],[tt';tt'],'Color',rgb(rgb22{i}),'LineWidth',2)
    %plot(p2*ones(size(idx2(tt)))-0.02,tt','.','color',[[218/255,165/255,32/255]]);
    plot(p2*ones(size(idx2(tt)))-0.02,tt','.','color',[[0.5,0.5,0.5]]);	
    max1 = max(max1,max(hyge_ind2(idx2(tt))));
    end
end

max2 = 0;
for i=1:22
    tt = find(chromoID_ind1(idx1)==i);
    if (nnz(tt)~=0)
    line([p1*ones(size(idx1(tt)));-(hyge_ind1(idx1(tt)))/2+p1],bb+[tt';tt'],'Color',rgb(rgb22{i}),'LineWidth',2)
    %plot(p1*ones(size(idx1(tt)))+0.02,bb+tt','.','color',[[218/255,165/255,32/255]]);
    plot(p1*ones(size(idx1(tt)))+0.02,bb+tt','.','color',[[0.5,0.5,0.5]]);	
    max2 = max(max2,max(hyge_ind1(idx1(tt)))); 
    end
end

[ii jj] = find(MM==1);
if riskornot==0
    %line([p1*ones(size(ii))';p2*ones(size(jj))'],[ii'+bb;jj'],'Color',rgb('LightCoral'),'LineStyle','-');
    %line([p1*ones(size(ii))'+0.04;p2*ones(size(jj))'-0.04],[ii'+bb;jj'],'Color',[173/255,216/255,230/255],'LineStyle','-');
     line([p1*ones(size(ii))'+0.04;p2*ones(size(jj))'-0.04],[ii'+bb;jj'],'Color',[0.82,0.82,0.82],'LineStyle','-','LineWidth',0.25);
else
    %line([p1*ones(size(ii))';p2*ones(size(jj))'],[ii'+bb;jj'],'Color',rgb('LightCoral'),'LineStyle','-');
    %line([p1*ones(size(ii))'+0.04;p2*ones(size(jj))'-0.04],[ii'+bb;jj'],'Color',[255/255,182/255,193/255],'LineStyle','-');
     line([p1*ones(size(ii))'+0.04;p2*ones(size(jj))'-0.04],[ii'+bb;jj'],'Color',[0.75,0.75,0.75],'LineStyle','-','LineWidth',0.25);	
end

s1 = p2+(-log10(0.05))/2;
s2 = -(-log10(0.05))/2;
line([s2 s2],[0 length(idx2)+1],'LineStyle','--','Color', [0.5 0.5 0.5])
line([s1 s1],[0 length(idx2)+1],'LIneStyle','--','COlor',[0.5,0.5,0.5])

%end
%xlim([-1 5])
ylim([-5 length(idx2)+5]);
% xlim([-1.5 3.5])
xlim([-max(hyge_ind1(ii))/2-0.1, max(hyge_ind2(jj))/2+2+0.1])
ax1=gca;

setfig
set(ax1,'box','off')
set(ax1,'Ycolor','w')
set(ax1,'tickdir','out')
% set(ax1,'Xtick',[-1.5:0.5:3.5],'XtickLabel',{'2.5','2','1','0','','','','0','1','2','2.5'})
xt = [round(-max(hyge_ind1(ii)))*100/100+1:0.5:round(max(hyge_ind2(jj)))*100/100+1];
xl = arrayfun(@(x)num2str(2*abs(x)),round(-max(hyge_ind1(ii)))*100/100:0.5:round(max(hyge_ind2(jj)))*100/100,'Uniform',0);
set(ax1,'Xtick',xt,'XtickLabel',xl)

set(gca,'YTick',[])
set(gca,'YColor','w')

t = title(sprintf('%s \n %s',pathname1,pathname2)) 
set(t,'Interpreter','none')

saveas(gcf,outputname)
