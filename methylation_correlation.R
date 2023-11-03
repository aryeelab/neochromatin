#################################################################
## Plot correlations between expression, methylation and coverage
library(rtracklayer)
library(GenomicRanges)

load('../Expression_data/DESeq_results.RData')
load('../Expression_data/lps853_DESeq_results.RData')

lps141.exp <- res2vs0$baseMean
names(lps141.exp) <- rownames(res2vs0)
lps853.exp <- lps853.res2vs0$baseMean
names(lps853.exp) <- rownames(lps853.res2vs0)

lps141.bg <- import('../Bedgraphs/LPS141_barcode01_5mc_m_combined.bedgraph')
lps853.bg <- import('../Bedgraphs/LPS853_barcode02_5mc_m_combined.bedgraph')


genes.gtf <- read.table('../hg38.refGene.gtf', sep='\t')
genes.gtf <- genes.gtf[genes.gtf$V3=='transcript',]
gene.names <- unlist(lapply(strsplit(genes.gtf$V9, split = ';'), function(x){gsub("gene_id ", "",  x[[1]])}))

genes.gtf <- as.data.frame(cbind(genes.gtf$V1, genes.gtf$V4, genes.gtf$V5, genes.gtf$V6, genes.gtf$V7, gene.names))

lps141.matrix <- matrix(ncol=3, nrow=length(genes.gtf$gene.names))
rownames(lps141.matrix) <- genes.gtf$gene.names
colnames(lps141.matrix) <- c('expression','methylation','coverage')

lps853.matrix <- matrix(ncol=3, nrow=length(genes.gtf$gene.names))
rownames(lps853.matrix) <- genes.gtf$gene.names
colnames(lps853.matrix) <- c('expression','methylation','coverage')


## This takes forever. Should have done full gr overlap and then subset
#  but hey, it works

for (gene in genes.gtf$gene.names){
  gene.gtf <- genes.gtf[genes.gtf$gene.names==gene,]
  gene.gr <- GRanges(seqnames = gene.gtf$V1, 
                     ranges = IRanges(as.numeric(gene.gtf$V2), as.numeric(gene.gtf$V3)))
  
  ## LPS141
  gene.meth.subset.gr <- subsetByOverlaps(lps141.bg, gene.gr)
  
  lps141.matrix[gene, 1] <- lps141.exp[gene]
  lps141.matrix[gene, 2] <- mean(gene.meth.subset.gr$score)
  lps141.matrix[gene, 3] <- mean(gene.meth.subset.gr$NA.)
  
  ## LPS853
  gene.meth.subset.gr <- subsetByOverlaps(lps853.bg, gene.gr)
  
  lps853.matrix[gene, 1] <- lps853.exp[gene]
  lps853.matrix[gene, 2] <- mean(gene.meth.subset.gr$score)
  lps853.matrix[gene, 3] <- mean(gene.meth.subset.gr$NA.)
}

save(lps141.matrix, lps853.matrix, file = '../R_objects/methylation_matrices.RData')


lps141.matrix.indices <- which(apply(lps141.matrix, 1, function(x){!any(is.na(x))}))
lps853.matrix.indices <- which(apply(lps853.matrix, 1, function(x){!any(is.na(x))}))

common.indices <- intersect(lps141.matrix.indices, lps853.matrix.indices)

lps141.common <- lps141.matrix[common.indices,]
lps853.common <- lps853.matrix[common.indices,]


rna.delta <- log(lps141.common[,1])-log(lps853.common[,1])
methylation.delta <- lps141.common[,2]-lps853.common[,2]
coverage.delta <- lps141.common[,3]-lps853.common[,3]


save(rna.delta, methylation.delta, coverage.delta, file = '../R_objects/delta_objects.RData')

library(scales)
pdf(file = '../plots/deltas_histograms.pdf', width = 5, height = 4)
par(mfrow=c(1,3))
hist(rna.delta, breaks=100, main='delta log rna'); hist(methylation.delta, breaks=100, main='delta meth'); hist(coverage.delta, breaks=100, main='delta coverage')
dev.off()

pdf(file = '../plots/deltas_correlations.pdf', width = 5, height = 4)
par(mfrow=c(1,3))
plot(coverage.delta, methylation.delta, frame.plot=F, pch=19, col=alpha('black', 0.3), main='delta coverage vs delta methylation')
plot(coverage.delta, rna.delta, frame.plot=F, pch=19, col=alpha('black', 0.3), main='delta coverage vs delta expression')
plot(methylation.delta, rna.delta, frame.plot=F, pch=19, col=alpha('black', 0.3), main='delta methylation vs delta expression')
dev.off()