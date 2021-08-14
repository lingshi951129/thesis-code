a<-read.table(file="E:/project/TCGA/same_result_500kb_wt_TPM_flt.txt",header=T,sep= "\t")
b<-read.table(file="E:/project/TCGA/dif_result_500kb_wt_TPM_flt.txt",header=T,sep= "\t")
a1<-a[,5]/1000
a2<-a[,6]
a3<-cbind(a1,a2)
a4<-a3[order(a3[,1]),]
aa<-data.frame(aa1=a4[,1],aa2=a4[,2])
aa$index <- 1:nrow(aa)
b1<-b[,5]/1000
b2<-b[,6]
b3<-cbind(b1,b2)
b4<-b3[order(b3[,1]),]
bb<-data.frame(bb1=b4[,1],bb2=b4[,2])
bb$index <- 1:nrow(bb)
#jpeg("fig1e_new.jpg")
loessMod10a <- loess(aa2 ~ index, data=aa, span=0.10)
smoothed10a <- predict(loessMod10a)
loessMod10b <- loess(bb2 ~ index, data=bb, span=0.10)
smoothed10b <- predict(loessMod10b)
plot(smoothed10a, x=aa$aa1, type="l",lwd=8,col="blue",xlab="Distance(kb)",ylab="Pearson Correlation",ylim=c(0,0.7),cex.axis=1.2,cex.lab=1.5)
lines(smoothed10b,x=bb$bb1,col="orange",lwd=8)
par(cex = 1.2)
legend("topright",c("Same domain","Cross boundary"),lty=1,col=c("blue","orange"))
