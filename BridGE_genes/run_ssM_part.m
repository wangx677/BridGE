function run_ssM_part(n1,n2)

for R=n1:n2
     computessi('RR',0,0.05,0.05,'plinkFile.cluster2',1,R);
     computessi('DD',0,0.05,0.05,'plinkFile.cluster2',1,R);
     computessi('RD',0,0.05,0.05,'plinkFile.cluster2',1,R);
end

for R=n1:n2
     computessi('combined',0,0.05,0.05,'plinkFile.cluster2',1,R);
end

