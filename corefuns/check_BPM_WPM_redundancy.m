function [BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group,WPM_group,PATH_group] = check_BPM_WPM_redundancy(fdrBPM,fdrWPM,fdrPATH,bpmindfile,FDRcut)

load(bpmindfile)
ssM{1} = 0;ssM{2}=0; % arbitrary -- it's not actully used because the overlap is based on pathway genes not actual interactions.  
j = 1;
for fdrcut=0.05:0.05:FDRcut
        ind = find(fdrBPM<=fdrcut);
        ind1 = ind(1:nnz(ind<=length(BPM.size)));
        ind2 = ind(nnz(ind<=length(BPM.size))+1:end)-length(BPM.size);

        ind = find(fdrWPM<=fdrcut);
	   ind3 = ind(1:nnz(ind<=length(WPM.size)));
        ind4 = ind(nnz(ind<=length(WPM.size))+1:end)-length(WPM.size);

        if exist('ind1')==0
                ind1=[];
        end

        if exist('ind2')==0
                ind2=[];
        end

        if exist('ind3')==0
                ind3=[];
        end

        if exist('ind4')==0
                ind4=[];
        end

        if (length(ind1)>1)
                BPM_sim{1} = bpmsim(BPM.ind1(ind1),BPM.ind2(ind1),BPM.ind1(ind1),BPM.ind2(ind1),2,ssM{1});
                TTT=BPM_sim{1}>=0.25;
                [noRD{1},group{1}] = graphconncomp(sparse(TTT));
                %[num2cell(group') path1(ind1)' path2(ind1)']
                a0 = length(unique(group{1}));
        elseif (length(ind1)==1)
                a0=1;
		group{1} = 1;
        else
                a0 = 0;
		group{1} = [];
        end

        if (length(ind2)>1)
                BPM_sim{2} = bpmsim(BPM.ind1(ind2),BPM.ind2(ind2),BPM.ind1(ind2),BPM.ind2(ind2),2,ssM{1});
                TTT=BPM_sim{2}>=0.25;
                [noRD{2},group{2}] = graphconncomp(sparse(TTT));
                %[num2cell(group') path1(ind2)' path2(ind2)']
                b0 = length(unique(group{2}));
        elseif(length(ind2)==1)
                b0 = 1;
		group{2} = 1;
        else
                b0 = 0;
		group{2} = [];
        end

        if (length(ind3)>1)
                BPM_sim{3} = bpmsim(WPM.ind(ind3),WPM.ind(ind3),WPM.ind(ind3),WPM.ind(ind3),2,ssM{1});
                TTT=BPM_sim{3}>=0.25;
                [noRD{3},group{3}] = graphconncomp(sparse(TTT));
                %[num2cell(group') path1(ind1)' path2(ind1)']
                c0 = length(unique(group{3}));
        elseif (length(ind3)==1)
                c0 = 1;
		group{3} = 1;
        else
                c0 = 0;
		group{3} = [];
        end

         if (length(ind4)>1)
                BPM_sim{4} = bpmsim(WPM.ind(ind4),WPM.ind(ind4),WPM.ind(ind4),WPM.ind(ind4),2,ssM{1});
                TTT=BPM_sim{4}>=0.25;
                [noRD{4},group{4}] = graphconncomp(sparse(TTT));
                %[num2cell(group') path1(ind1)' path2(ind1)']
                d0 = length(unique(group{4}));
        elseif (length(ind4)==1)
                d0 = 1;
		group{4} = 1;
        else
                d0 = 0;
		group{4} = [];
        end

	g1 = length(unique(group{1}));
	g2 = length(unique(group{2}));
	g3 = length(unique(group{3}));
	g4 = length(unique(group{4}));

	BPM_group{j} = [group{1} group{2}+g1];
	BPM_nosig_noRD(j) = a0 + b0;
	
	WPM_group{j} = [group{3} group{4}+g3];
	WPM_nosig_noRD(j) = c0 + d0;

	clear ind1 ind2 ind3 ind4 ind noRD group TTT BPM_sim a0 b0 c0 d0	
	% Pathways
	ind = find(fdrPATH<=fdrcut);
	ind1 = ind(1:nnz(ind<=length(WPM.size)));
	ind2 = ind(nnz(ind<=length(WPM.size))+1:end)-length(WPM.size);

	if exist('ind1')==0
                ind1=[];
        end

        if exist('ind2')==0
                ind2=[];
        end

	if (length(ind1)>1)
		PATH_sim{1} = pathsim(WPM.ind(ind1));
		TTT = PATH_sim{1}>=0.25;
		[noRD{1},group{1}] = graphconncomp(sparse(TTT));
		a0 = length(unique(group{1}));
	elseif (length(ind1)==1)
		a0 = 1;
		group{1} = 1;
	else
		a0 = 0;
		group{1} = [];
	end


	if (length(ind2)>1)
                PATH_sim{2} = pathsim(WPM.ind(ind2));
                TTT = PATH_sim{2}>=0.25;
                [noRD{2},group{2}] = graphconncomp(sparse(TTT));
                b0 = length(unique(group{2}));
        elseif (length(ind2)==1)
                b0 = 1;
		group{2} = 1;
        else
                b0 = 0;
		group{2} = [];
        end
	
	g1 = length(unique(group{1}));
        g2 = length(unique(group{2}));

	PATH_group{j} = [group{1} group{2}+g1];
	PATH_nosig_noRD(j) = a0 + b0;
	j = j+1;
end

