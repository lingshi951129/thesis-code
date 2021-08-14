library(statmod)
library(edgeR)
a <- read.csv(file='E:/wt.raw.txt', header= T, sep='\t')
b <- read.csv(file='E:/mut.raw.txt', header= T, sep='\t')
#c <- read.csv(file='E:/project/K562/RNA_Seq/RHC9922.htseq.counts', header= F, sep='\t')
#d <- read.csv(file='E:/project/K562/RNA_Seq/RHC9923.htseq.counts', header= F, sep='\t')
#data <- cbind(a[,-1],b[,-1],c[,-1],d[,-1])
data <- cbind(a[,-1],b[,-1])
rownames(data) <- a[,1]
#condi <- c('Con1','Con2','KO_1','KO_2')
#colnames(data) <- condi
condi<-colnames(data) 
group <- factor(substring(condi,1,3))
group <- relevel(group, ref="wt_")
time <- factor(substring(condi,4,6))
time
y <- DGEList(counts=data,group=group)
y
keep <- rowSums(cpm(y)>1) >= 2
keepexp <- filterByExpr(y, group=group)
y <- y[keepexp, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~time+group)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- topTags(glmQLFTest(fit), n = nrow(y$counts))
write.table(qlf, 'E:/DNMT3A_raw_DEG1.txt', sep = '\t', col.names = NA, quote = FALSE)

fig_data <- read.delim('E:/DNMT3A_raw_DEG1.txt',row.names = 1, sep='\t',check.names = F) 
fig_data <- fig_data[order(fig_data$FDR,fig_data$logFC,decreasing = c(F,T)),]
fig_data[which(fig_data$logFC >= 1 & fig_data$FDR < 0.05),'sig'] <- 'up'
fig_data[which(fig_data$logFC <= -1 & fig_data$FDR < 0.05),'sig'] <- 'down'
fig_data[which(abs(fig_data$logFC) <= 1 | fig_data$FDR >= 0.05),'sig'] <- 'none'
library(ggplot2)
library(extrafont)
loadfonts(device = "win")
p <- ggplot(data=fig_data, aes(x=logFC, y=-log10(FDR),color = sig))+
  geom_point(size=1)+
  scale_color_manual(values=c('orangered2','gray','palegreen3'),limits=c('up','none','down'))+
  labs(x='log2 Fold Change',y='-log10 FDR',title='DNMT3A mutant vs Wild Type',color='')+
  theme(text=element_text(size=16,family="Arial"),plot.title=element_text(hjust=0.5,size=18),panel.grid=element_blank(),
        panel.background=element_rect(color='black',fill='transparent'),
        legend.key=element_rect(fill='transparent'))+
  geom_vline(xintercept=c(-1,1),lty=3,color='black')+
  geom_hline(yintercept=2,lty=3,color='black')+
  xlim(-20,20)+ylim(0,35)
p