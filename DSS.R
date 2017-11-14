.libPaths( c( .libPaths(), "/home/meng_yang/R/x86_64-pc-linux-gnu-library/3.3"))
library("DSS")
library("reshape2")
args <- commandArgs(TRUE)
#### usage: Rscript script.R  comparegroup outname
data_path = args[1] # count
compareGroup = args[2] # LC002-093,LC002-094--LC002-095,LC002-096
compareGroupname = args[3] # aa--bb
out=args[4]
print(compareGroup)

groupname = unlist(strsplit(compareGroupname,"--"))
sampleName  = unlist(strsplit(compareGroup,"--|,"))
compareGroupNames = unlist(strsplit(compareGroup,"--"))
sampleName1 = unlist(strsplit(compareGroupNames[1],","))
sampleName2 = unlist(strsplit(compareGroupNames[2],","))

outdir = paste(out,compareGroupname,sep="/")
dir.create(outdir,recursive = TRUE)
outname = paste(outdir,compareGroupname,sep="/")
x = lapply(sampleName, function(f){read.table(list.files(path=data_path, pattern = f, full.names = TRUE), header = TRUE,check.names = FALSE)})
BSobj <-  makeBSseqData(x,sampleName)
#========================== output methdata matrix ======================#
N_m = getBSseq(BSobj,type=c("Cov"))
X_m = getBSseq(BSobj,type=c("M"))
chr_m = as.data.frame(getBSseq(BSobj,type=c("gr")))[,1:2]
N_m.melt = melt(
  cbind(pos = paste(chr_m$seqnames, chr_m$start, sep=":"),
        as.data.frame(N_m))
         , id.vars = 1, measure.vars=2:(length(sampleName)+1),value.name = "coverage", variable.name = "sample")
X_m.melt = melt(
  cbind(pos = paste(chr_m$seqnames, chr_m$start, sep=":"),
        as.data.frame(X_m))
        , id.vars = 1, measure.vars=2:(length(sampleName)+1),value.name = "methylated", variable.name = "sample")
tmp.data <- merge(N_m.melt, X_m.melt, by = c("pos","sample"))
tmp.melt <- melt(tmp.data, id.vars = 1:2, measure.vars = 3:4, variable.name = "variable", value.name = "value")
tmp <- dcast(tmp.melt, pos ~ sample + variable)
loc = data.frame(pos=paste(chr_m$seqnames, chr_m$start,sep=":"))
meth.data=tmp[loc$pos,]
write.table(meth.data, file=paste(outname,".methdata.matrix.xls",sep=""),sep="\t",row.names=F,quote=F)

methylevel=X_m/N_m
#methylevel.data=data.frame(chr_m,methylevel)
methylevel.data=data.frame(loc,methylevel,check.names=FALSE)
write.table(methylevel.data, file=paste(outname,".methlevel.matrix.xls",sep=""),sep="\t",row.names=F,quote=F)

#========= filter out CpGs of low coverage & few samples =================#
X = getBSseq(BSobj, "M")
N = getBSseq(BSobj, "Cov")
N[N<10]  = 0
X[N<10]  = 0
select.loc = which((rowSums(N[,sampleName1] >= 10) >= 0.8*length(sampleName1)) & (rowSums(N[,sampleName2] >= 10) >= 0.8*length(sampleName2)))
BSobj=BSobj[select.loc,]
print(length(N[,1]))
print (length(select.loc))
X = getBSseq(BSobj, "M")
N = getBSseq(BSobj, "Cov")
chr = as.data.frame(getBSseq(BSobj,type=c("gr")))[,1:2]
meancov=data.frame(pos=paste(chr$seqnames, chr$start,sep=":"),meanmethcov1=rowMeans(N[,sampleName1], na.rm=TRUE),meanmethcov2=rowMeans(N[,sampleName2], na.rm=TRUE))
meancov$meanmethcov1=round(meancov$meanmethcov1)
meancov$meanmethcov2=round(meancov$meanmethcov2)

#=========DMLtest=========================================================#
dmlTest <- DMLtest(BSobj, group1=sampleName1, group2=sampleName2)
dmlTest.pos = data.frame(pos=paste(dmlTest$chr,dmlTest$pos,sep=":"),dmlTest)
dmlTest.out = merge(meancov,dmlTest.pos,by="pos")
dmlTest.out = dmlTest.out[,c(-4,-5)]
colnames(dmlTest.out) = c("chr:pos",paste(groupname[1],".meancov",sep=""),paste(groupname[2],".meancov",sep=""),paste(groupname[1],".mu",sep=""),paste(groupname[2],".mu",sep=""),"diff","diff.se","stat",paste(groupname[1],".phi",sep=""),paste(groupname[2],".phi",sep=""),"pval","fdr")
write.table(dmlTest.out, file=paste(outname,".DMLtest.xls",sep=""),sep="\t",row.names=F,quote=F)

#=============diff & p plot===============================================#
mygrey <- rgb(89/255,87/255,87/255)
#volcano plot
png(file=paste(outname,".volcano_plot.png",sep=""),width = 480, height = 480)
par(family="sans",col.main=mygrey,col=mygrey,col.lab=mygrey,col.axis=mygrey,cex.main=2,cex.lab=1.5,cex.axis=1.5,mgp=c(2.5,1,0))
dmlTest$logp = -log10(dmlTest$pval)
dmlTest$logfdr = -log10(dmlTest$fdr)
plot(dmlTest[,"diff"],dmlTest[,"logp"],mgp=c(2.5,1,0),xlab="diff",ylab="-log10(pvalue)",cex.lab =1.7,pch=20,cex=0.5,col=rgb(197/255,27/255,138/255),main="Volcano Plot",xlim=range(dmlTest[,"diff"]))
box()
dev.off()
#p-value histogram
png(file=paste(outname,".p-value.png",sep=""),width = 480, height = 480)
par(family="sans",col.main=mygrey,col=mygrey,col.lab=mygrey,col.axis=mygrey,cex.main=2,cex.lab=1.5,cex.axis=1.5,mgp=c(2.5,1,0))
hist(dmlTest$pval,col=rgb(250/255,159/255,181/255),xlab="p-value",ylab="Frequency",main="Histogram Of p-value")
box()
dev.off()

#======================fun======================
DMR_MethylRatio <- function (OneDMR, BSobj, ext = 500, ylim = c(0,1), plot=FALSE){
  allchr = as.character(seqnames(BSobj))
  allpos = start(BSobj)
  X = getBSseq(BSobj, "M")
  N = getBSseq(BSobj, "Cov")
  chr = as.character(OneDMR$chr)
  ix.chr = which(allchr == chr)
  thispos = allpos[ix.chr]
 thisN = N[ix.chr, ]
  thisX = X[ix.chr, ]
 xlim = c(OneDMR$start - ext, OneDMR$end + ext)
 ix1 = which(thispos <= xlim[2] & thispos >= xlim[1])
  ix2 = which(thispos <= OneDMR$end  & thispos >= OneDMR$start )

  thisP = thisX/thisN
  if (plot){
    nSample = ncol(X)
   if (nSample > 2) {
     y.cex = 0.66
    }
    else y.cex = 1
   sNames = sampleNames(BSobj)
    par(mfrow = c(nSample, 1), mar = c(2.5, 2.5, 1.6, 2.5), mgp = c(1.5,0.5, 0))
    for (i in 1:ncol(X)) {
      plot(thispos[ix1], thisP[ix1, i], type = "h", col = "blue",
          axes = F, lwd = 1.5, xlab = "", ylab = "", ylim = ylim,
           xlim = xlim, main = sNames[i])
      box(col = "black")
      axis(1, )
      axis(2, col = "blue", col.axis = "blue")
      mtext(chr, side = 1, line = 1.33, cex = y.cex)
      mtext("methyl%", side = 2, line = 1.33, col = "blue",
            cex = y.cex)
      thisN.norm = thisN[ix1, i]/max(thisN[ix1, ]) * ylim[2]
      lines(thispos[ix1], thisN.norm, type = "l", col = "gray",
            lwd = 1.5)
      axis(side = 4, at = seq(0, ylim[2], length.out = 5),
           labels = round(seq(0, max(thisN[ix1, ]), length.out = 5)))
      mtext("read depth", side = 4, line = 1.33, cex = y.cex)
     rect(OneDMR$start, ylim[1], OneDMR$end, ylim[2], col = "#FF00001A",
           border = NA)
    }
    return(thisP[ix2, ])
  } else {
    return(colMeans(thisP[ix2, ], na.rm=TRUE))
  }
}

#====================callDML of delta & p =============================

#for (i in c(0.6,0.7,0.75,0.8,0.85,0.9)){
for (i in c(0.009,0.01,0.02,0.03,0.04,0.05)){
	for(j in c(0.2,0.1,0.05)){
#	cutoff_i=quantile(abs(data.frame(dmlTest)[,"diff"]),i)
#	cutoff_i=round(cutoff_i,digits=4)
	cutoff_i = i
	tmp_i=i
	tmp_i=paste(i,cutoff_i,sep="-")
	tmp = paste(tmp_i,j,sep="-");
	tmp0 = paste(".dml.",tmp,sep="")
	tmp = paste(".dmr.",tmp,sep="")
	tmp0 = paste(tmp0,".xls",sep="")
	tmp = paste(tmp,".xls",sep="")
	dmls <- callDML(dmlTest,delta=cutoff_i,p.threshold=j)
	dmrs <- callDMR(dmlTest,delta=cutoff_i,p.threshold=j)
	flag = try(y <- lapply(
	  split(dmrs, paste(dmrs$chr, dmrs$start, sep = ":")),
	  function(dm){DMR_MethylRatio(dm,BSobj)}
	),TRUE)
	if (class(flag)=="try-error"){next}
	dmrs_new  <- cbind(
		dmrs,
	  	pos = paste(dmrs$chr, dmrs$start, sep = ":")
	)
	for (k in 1:length(y[[1]])){
		dmrs_new <- cbind(
		dmrs_new,
		V1 = unlist(lapply(y,function(x){ as.numeric(x[k])}))[dmrs_new$pos]
	  )
	  names(dmrs_new) <- c(  names(dmrs_new)[-c(dim(dmrs_new)[2])], names(y[[1]])[k] )
	}
	dmls.pos = data.frame(pos=paste(dmls$chr,dmls$pos,sep=":"),dmls)
	dmls.out = merge(meancov,dmls.pos,by="pos")
	dmls.out = dmls.out[,c(-4,-5,-15,-16)]
	colnames(dmls.out) = c("chr:pos",paste(groupname[1],".meancov",sep=""),paste(groupname[2],".meancov",sep=""),paste(groupname[1],".mu",sep=""),paste(groupname[2],".mu",sep=""),"diff","diff.se","stat",paste(groupname[1],".phi",sep=""),paste(groupname[2],".phi",sep=""),"pval","fdr","postprob.overThreshold")
	
	colnames(dmrs_new) = sub("meanMethy1",paste(groupname[1],"meanMethy",sep="."),colnames(dmrs_new))
	colnames(dmrs_new) = sub("meanMethy2",paste(groupname[2],"meanMethy",sep="."),colnames(dmrs_new))
	write.table( dmls.out, file=paste(outname,tmp0,sep=""),sep="\t",row.names=F,quote=F)
	write.table( dmrs_new, file=paste(outname,tmp,sep=""),sep="\t",row.names=F,quote=F)
	}
}



