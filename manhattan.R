chrs <- snakemake@params[["chromosomes"]]
##cat ("chromosomes are:", chrs, "\n")
cat ("inputfiles are \n")
infiles=snakemake@input[["infiles"]]
for (file in infiles){
    cat (file, "\n")
}

plot_file=snakemake@output[["plot_file"]]
#txt_file=snakemake@output[["txt_file"]]


cat ("parameters\n")
##parameters
#r2_thresh <- as.numeric(snakemake@params[["r2_thresh"]])
maf_thresh <- as.numeric(snakemake@params[["maf_thresh"]])

cat ("maf thresh is", maf_thresh, "\n")
#cat ("r2 thresh is", r2_thresh, "\n")



##info
cat ("infos\n")
phenotype=snakemake@wildcards[["cat"]]
sex=snakemake@wildcards[["sex"]]

mytitle=paste0(phenotype, "_", sex)
cat ("mytitle is ", mytitle, "\n")

cnames  <- c('chr', 'pos', 'af', 'p_wald', 'p_lrt', 'p_score')


for (i in 1:length (infiles)) {
    inf  <- matrix(   scan (infiles [i]),ncol=6,byrow=T )
    maf  <- pmin (inf [,3], 1-inf[,3])
    finf  <- inf [maf >= maf_thresh,]

    ##add chr column
    ##finf  <- cbind (rep(chrs[i], nrow(finf)), finf)
    if (i ==1) {
        allchr <- finf
    }else{
        allchr  <- rbind (allchr, finf)    
    }
    
}


res  <- data.frame (allchr)
head (res)


tiff (plot_file, height=12,width=20,units="in",res=300)
chrcol <- 1
poscol <- 2
pcol <- 4
colors <- c("steelblue", "indianred2")
fontsize=1.5  #font size for the ylab
#leave larger margins for ylab
mymar <- c (4, 6, 2, 2)
par (mar=mymar, cex.lab=1.5,cex.axis=1.5)
thresh <- 6
#-------------------------------------------------

#alternate the color for consecutive chromosomes
mycol = rep (NA, nrow (res))
evens <- chrs [1:length(chrs)%%2==0]
odds <- chrs [1:length(chrs)%%2==1]
mycol [which (res[,chrcol] %in% evens) ]  <-  colors [1]
mycol [which (res[,chrcol] %in% odds) ]  <-  colors [2]

##mycol [  res [,3] <= maf_thresh | res [,4] <= r2_thresh ]  <- "white"

#get the chromosome sizes
##chrs <- sort(unique (res [, chrcol]))
sizes <- rep(0,29)
for (i in 1:length (chrs)) {
    mychr  <- res [res [,chrcol] == chrs[i], ]
    if (nrow (mychr) > 0) {
        sizes [chrs [i]]  <- max (mychr[,poscol])
    }else {
        sizes [chrs [i]]  <-  0
    }
}

#continuous positions
toadd <- c (0, cumsum (sizes))
toadd <- toadd [-c(length(toadd))]
res$pos <- res[, poscol] + toadd [res [, chrcol]]


cat ("myrange :", min (res[,pcol]), max(res[,pcol]), "\n" )
maxi <- max (-log10(res[, pcol]))
plot (res$pos, -log10(res[, pcol]), col=mycol, xaxt="n", ylab=expression (-log [10](italic("P"))), ylim =c(0, maxi*1.1), cex.main=1 , cex.lab=fontsize, cex.axis=1, xlab="", pch=20, cex=1.5,yaxt="n")
axis (side=2, las=2)

#position for chr label.. midpoint
##ats <- (sizes/2) + toadd
##oats  <- ats
ats  <- c ()
for (chr in chrs) {
   myats  <-  mean (range (res [res[,1] == chr,"pos"]))
   ats  <- c (ats, myats)
   cat (chr, "\n")
}



#axis usually no space for axis
#axis (side=1,  at=ats, labels=1:29)
text (ats, rep(maxi*1.1, length(ats)), labels=chrs, cex=1)

#vertical lines separting the chromosomes
abline (v= sizes + toadd,col='gray', lty=2 )

#threshold for significance
abline (h =thresh,  col='darkseagreen', lty=2 )

dev.off ()

