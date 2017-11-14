args<-commandArgs(TRUE)
library(ggplot2)
input=args[1]
out=args[2]
flag=args[3]

data <- read.table(input,header=T,sep="\t",check.names=FALSE)
newdata <- subset(data,!grepl("all",anno.gene.type))
mylabel <- paste(gsub(" ","",newdata$anno.gene.type),"(",gsub(" ","",newdata$percent),")")
newdata$fill = factor(rownames(newdata),levels=rownames(newdata))
mylabel <- paste(gsub(" ","",newdata$anno.gene.type),"(",gsub(" ","",newdata$percent),")")
colours <- c(rgb(192/255,27/255,125/255),rgb(233/255,163/255,201/255),rgb(253/255,224/255,239/255),rgb(230/255,245/255,208/255),rgb(161/255,215/255,106/255),rgb(77/255,146/255,33/255), gray.colors(10))
flag=paste("                        ",paste(flag,"Annotation",sep=" "),sep=" ")
print(flag)
pdf(file=out, width=10, height=10)
p <- ggplot(newdata,aes(x="",y=count,fill=fill))
p = p + geom_bar(stat = "identity",width=1,colour="white") 
p = p + coord_polar(theta="y") 
p = p + scale_fill_manual(breaks=newdata$fill,values=colours,labels=mylabel) 
p = p + labs(x="", y="",title=flag) 
p = p + theme(axis.ticks=element_blank(),legend.title=element_blank(),axis.text.x=element_blank(),plot.title=element_text(size=20,face="bold",family="sans",colour=rgb(89/255,87/255,87/255)),panel.background=element_blank(), legend.text=element_text(size=15,family="sans",colour=rgb(89/255,87/255,87/255)))
p
dev.off()
