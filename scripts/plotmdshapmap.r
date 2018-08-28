#!/soft/r/3.2.2/linux_x86_64/bin/Rscript 

args <- commandArgs(trailingOnly = TRUE)
print(args)

mdsfile <- args[1]
popfile <- args[2]
plotfile <- args[3]
pos <- args[4]
pop <- args[5]
hapmappop <- args[6]
rm(args)

print(hapmappop)

if (is.na(hapmappop)) {
	hapmappop<-"GIH"
	}

print(hapmappop)

mds <- read.table(sprintf("%s.mds",mdsfile), header=T) 
popnids <- read.table(sprintf("%s",popfile), header=T)
pmsize<-dim(popnids)
mm <- matrix(0,pmsize[1],pmsize[2]-2)

#colorlist<-c("#e31a1c","#b15928","#ffff99","#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
colorlist<-c("black","red","green","blue","yellow","#ff1493","#ae017e","#dd3497","#f768a1","#fa9fb5")

for (i in 3:pmsize[2]) {
	v <- i-2
	mm[,v] = v
	}

colourM <- popnids[,3:pmsize[2]]*mm
colour <- rowSums(colourM)

pdf(sprintf("%s.pdf",plotfile))
plot(mds$C1, mds$C2, col=colorlist[colour], pch="+", xlab="PC1", ylab="PC2")
if (pop=="CEU"){
	ind<-which(colour %in% 1)
	points(mds$C1[ind], mds$C2[ind], col = colorlist[colour[ind]], pch = "+")
	} else if (pop=="CHB"){
	ind<-which(colour %in% 2)
	points(mds$C1[ind], mds$C2[ind], col = colorlist[colour[ind]], pch = "+")
	} else if (pop=="ASW"){
	ind<-which(colour %in% 3)
        points(mds$C1[ind], mds$C2[ind], col = colorlist[colour[ind]], pch = "+")
	 } else if (pop==hapmappop){
	ind<-which(colour %in% 4)	
        points(mds$C1[ind], mds$C2[ind], col = colorlist[colour[ind]], pch = "+")
	 } else if (pop=="YRI"){
	ind<-which(colour %in% 5)
        points(mds$C1[ind], mds$C2[ind], col = colorlist[colour[ind]], pch = "+")
	}

ll<-as.matrix(read.table(sprintf("%s",popfile),nrow=1))
ll<-ll[3:pmsize[2]]
if (hapmappop=="GIH"){
	ll[1:5]<-c("European","Han Chinese","African","Indian","Yoruba")
	} else if (hapmappop=="JPT") {
	ll[1:5]<-c("European","Han Chinese","African","Japanese","Yoruba")
	 } else if (hapmappop=="MEX") {
	ll[1:5]<-c("European","Han Chinese","African","Mexican","Yoruba")
	}

tt<-pmsize[2]-2
c1<-rep(1,tt)
c2<-rep(2.5,tt)
c3<-colorlist[1:tt]
# legend(pos, c("CEU","CHB","ASW","GIH","YRI","StudyPop"), lty=c1, lwd=c2, col=c3)
# legend(pos, c("European","Han Chinese","African","Indian","Yoruba","StudyPop"), lty=c1, lwd=c2, col=c3)
legend(pos, ll, lty=c1, lwd=c2, col=c3)
title("MDS plot colored by populations")

dev.off()
