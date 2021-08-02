function run_ssM

for R=0:1000
     computessi('RR',0,0.05,0.05,'plinkFile.cluster2',1,R);
     computessi('DD',0,0.05,0.05,'plinkFile.cluster2',1,R);
     computessi('RD',0,0.05,0.05,'plinkFile.cluster2',1,R);
end

for R=0:1000
     computessi('combined',0,0.05,0.05,'plinkFile.cluster2',1,R);
end

