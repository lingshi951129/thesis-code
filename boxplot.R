library(ggplot2)
data<-read.csv("G:/DNMT3A/DNMT3A_wt_mut.txt",header=T,sep = "\t")
p<-ggplot(data=data, aes(x=Types,y=Value))+geom_boxplot(aes(fill=Types),width=0.2,outlier.colour="red")+geom_point(size=1)+ggtitle("TCGA_LAML samples")+ylab("DNMT3A TPM")+xlab("")+theme(legend.position="right")
#+geom_jitter(shape=16, position=position_jitter(0.2))
p+ theme_bw()+  theme(text=element_text(size=16,family="Arial",face="bold"),plot.title=element_text(hjust=0.5,size=18,family = "Arial",face = "bold"),panel.grid=element_blank(),
                      panel.background=element_rect(color='black',fill='transparent'),
                      legend.key=element_rect(fill='transparent'),axis.text.x = element_text(size = 15, family = "Arial",face = "bold.italic"))
wt<-data[which(data$Types=="DNMT3A_wt"),]
a<-boxplot.stats(wt$Value)$stat
a
mut<-data[which(data$Types=="DNMT3A_mut"),]
b<-boxplot.stats(mut$Value)$stat
b
#p+geom_hline(yintercept=b,linetype="dotted")
#+annotate("text", x = 0.5, y = b[1], label = b[1])+annotate("text", x = 0.5, y = b[2], label = b[2])+annotate("text", x = 0.5, y = b[3], label = b[3])+annotate("text", x = 0.5, y = b[4], label = b[4])+annotate("text", x = 0.5, y = b[5], label = b[5])