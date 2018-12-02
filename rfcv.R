library(randomForest)
a <-read.table("genus.select.sub.tran.xls",head=T,sep="\t")
a<-a[-1]
c<- c(1:19)
colnames(a)<-c
b<-read.table("map.LT_S_LT_L.txt",head=T,sep="\t",,comment.char = "")
result <- rfcv(a, b$group, cv.fold=10,step=0.99, mtry=function(p) c(1: floor(sqrt(p))))
pdf("LT_S_LT_L.RFDV.pdf")
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2,col="red"))
dev.off()
write.table(file="LT_S_LT_L.error.xls",result$error.cv,sep="\t")
