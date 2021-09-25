#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(csaw))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(DESeq2))

parser <- ArgumentParser(description=
  'call peaks using DESeq2 packages.')

parser$add_argument("--treat", required=T, nargs="+", help="Bamfiles of sample in Treat group.")
parser$add_argument("--control", required=T, nargs="+", help="Bamfiles of sample in Control group.")
parser$add_argument("--sites", required=T, help="Filtered candidate sites.")
parser$add_argument("--window", required=F, type="integer", default=30, help="In which size of region to count reads.")
parser$add_argument("--out", required=T, help="File name to save Rdata")

args <- parser$parse_args()

# parse arguments
tbams <- args$treat
cbams <- args$control

condition <- factor(c(rep("control", length(cbams)), rep("treat", length(tbams))), levels=c("control","treat"))
ext.info <- c(rep(300, length(cbams)), rep(NA, length(tbams)))

# prepare windows
window <- args$window
half.window <- round(window/2, 0)
sites.df <- read.table(args$sites, sep='\t',header=T)
sites.df$start <- sites.df$position - half.window
sites.df$end <- sites.df$position + half.window
region.info <- paste(sites.df$chr, sites.df$start, sites.df$end, sep = '_')

site.regions <- GRanges(seqnames = sites.df$chr,
  ranges = IRanges(start = sites.df$start, end = sites.df$end))

# 1 step.
## make readcounts matrix
rp <- readParam(pe='none', dedup=FALSE)

exp <- regionCounts(c(cbams, tbams),regions=site.regions, ext=list(ext.info, NA), param=rp)
norm.factors <- normFactors(exp, se.out=F)

mtx <- exp@assays@data@listData$counts
exprSet <- as.data.frame(mtx)
colnames(exprSet) <- exp$bam.files
rownames(exprSet) <- region.info

lib.info <- data.frame(lib.sizes=exp$totals,
  row.names=exp$bam.files,
  norm.factors=norm.factors)

#save(exprSet, condition, lib.info, file=args$out)
coldata <- data.frame(row.names = colnames(exprSet), condition)
dds <- DESeqDataSetFromMatrix(countData = exprSet, colData = coldata, design= ~ condition )
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition','treat','control'))
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize=TRUE)), by="row.names", sort=FALSE)

write.table(resdata, file=args$out, sep='\t', col.names=T, row.names=T, quote=F)

cat("\nDone.")
