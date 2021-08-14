library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/home/lingshi/genome/hg38_p2/gencode.v22.chr_patch_hapl_scaff.annotation.gtf",format="gtf")
exons_gene <- exonsBy(txdb,by="gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_lens1 <- as.data.frame(exons_gene_lens)
exons_gene_lens1 <- t(exons_gene_lens1)
counts = read.csv("/home/lingshi/AML/RNA-Seq/metadata/merge.htseq.counts.txt",header = T,sep = "\t")
colnames(exons_gene_lens1)="Length"
count_with_length <- cbind(counts,exons_gene_lens1)
count_with_length <- count_with_length[,-1]
kb <- count_with_length$Length/1000
countdata <- count_with_length[,1:9]
rpk <- countdata/kb
fpkm <- t(t(rpk)/colSums(countdata)*10^6)
write.table(fpkm,file="/home/lingshi/AML/RNA-Seq/metadata/merge.FPKM.counts",sep="\t",quote=F)
