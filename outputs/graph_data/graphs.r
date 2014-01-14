getwd()
#pathtogff <- paste(getwd(),"/gffs",sep="") 
#setwd(pathtogff)
setwd("C:/Users/Jia/Google Drive/Academic/dis/src/outputs/graph_data/gffs")

getwd()
source("http://chromatin.bio.fsu.edu/DennisLab.R")
data1<-read.table("../baseseq.gff")
extfiles<-list.files(pattern=".gff")

pdf("../output.pdf")
for(f in extfiles) {
	data2<-read.table(f)
	drawGffPlots(data1,data2,trackingDye=TRUE,legend=TRUE,name1="base sequence", name2=f,newDev=0)
}
dev.off()