same <- read.table(file = "E:/project/TCGA/without_flt/same_delta_500kb_TPM_wo_flt.txt",header = F,sep = "\t",stringsAsFactors = F)
dif <- read.table(file = "E:/project/TCGA/without_flt/dif_delta_500kb_TPM_wo_flt.txt",header = F,sep = "\t",stringsAsFactors = F)
samed <- same[,3]
samep <- same[,4]
difd <- dif[,3]
difp <- dif[,4]
#("/home/lingshi/correlation/metadata/Volcano_plot_filtered.jpg")
plot(x=samed, y=samep, log="y", xlim=c(-1,1),ylim=c(1,1e-16),xaxt="n",yaxt="n",cex.lab = 1.5,                                                                                                xlab="Delta Correlation(DNMT3A mutant-DNMT3A wt)",ylab="Significance",type="p",col="darkgray"                                                                                                ,pch=20, bty="l",family="Calibri Light")
pvec <- 10^seq(-15,0,by=5)
axis(side=2,at=pvec,labels=pvec)
points(x=difd, y=-difp,log="y",col="darkgrey",pch=20)
points(samed[samep<1e-9],samep[samep<1e-9],col="blue",pch=20)
points(difd[difp<1e-9],difp[difp<1e-9],col="orange",pch=20)
legend("bottomright",c("Same domain","Cross boundary"),pch=18,col=c("blue","orange"),box.lty=                                                                                                0,bg=NULL,cex=1)
axis(side=1)
