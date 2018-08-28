#!/soft/r/3.2.2/linux_x86_64/bin/Rscript 

args <- commandArgs(trailingOnly = TRUE)
print(args)

mdsfile <- args[1]
plotfile <- args[2]
pos <- args[3]
rm(args)

mds <- read.table(sprintf("%s.mds",mdsfile), header=T) 
fam <- read.table(sprintf("%s.fam",mdsfile), header=F, col.names = c("FID","IID","PID","MID","sex","pheno"))

pdf(sprintf("%s.pdf",plotfile))

if ( all(sort(unique(fam$pheno)) == c(1,2)) ) {
        plot(mds$C1, mds$C2, col = fam$pheno, pch = "+", xlab = "PC1", ylab = "PC2")
	legend(pos, c("Controls","Cases"), lty = c(1,1), lwd = c(2.5,2.5), col = c(1,2))
        title("MDS plot colored by phenotype")
    }

if ( all(sort(unique(fam$sex)) == c(1,2)) ) {
	plot(mds$C1, mds$C2, col = fam$sex, pch = "+", xlab = "PC1", ylab = "PC2")
	legend(pos, c("Male","Female"), lty = c(1,1), lwd = c(2.5,2.5), col = c(1,2))
	title("MDS plot colored by sex")
    }

dev.off()
