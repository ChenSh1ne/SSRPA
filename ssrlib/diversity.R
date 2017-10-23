Args <- commandArgs()
ssrfile <- Args[6]
ploidy <- Args[7]
outpath <- Args[8]

library(adegenet)
library(poppr)
library(mmod)
#library(ploysat)
library(pegas)
library(hierfstat)
library(phangorn)

# genetic diversity for each single loci
obj_genind <- read.genalex(ssrfile, ploidy = as.numeric(ploidy), geo =FALSE, region=FALSE,genclone=FALSE)
div <- summary(obj_genind)
div_frame <- data.frame(div$loc.n.all, div$Hobs, div$Hexp)
outfile = paste(outpath,'/diversity.txt',sep='')
write.table(div_frame, file=outfile)

summ_stats <- diff_stats(obj_genind)
outfile2 <- paste(outpath,'/summ_stats.locus',sep='')
outfile3 <- paste(outpath,'/summ_stats.global',sep='')
write.table(summ_stats$per.locus, file=outfile2)
write.table(summ_stats$global, file=outfile3)


# group differenciation
pop_poppr <- poppr(obj_genind)
outfile4 <- paste(outpath,'/poppr.txt',sep='')
write.table(pop_poppr, file=outfile4)

#PCA
X <- scaleGen(obj_genind, NA.method="mean")
pca <- dudi.pca(X, cent=FALSE, scale=FALSE, scannf=FALSE, nf=3)
outfile5 <- paste(outpath, '/pca.pdf',sep='')
pdf(file=outfile5)
col <- funky(15)
s.class(pca$li, pop(obj_genind), col=transp(col, .6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
title("PCA Analysis")
dev.off()

# pairwise Fst, Gst...
matFST <- pairwise.fst(obj_genind, res.type="matrix")
outfile6 <- paste(outpath, '/fst.txt', sep='')
write.table(matFST, file=outfile6)

# pop tree
popobj <- genind2genpop(obj_genind)
popdist <- dist.genpop(popobj)
poptree <- upgma(popdist)
outfile7 <- paste(outpath, '/poptree.pdf',sep='')
pdf(file=outfile7)
plot(poptree)
dev.off()

#dendrogram

