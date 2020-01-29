#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#print(args[4])
print(getwd())
print(list.files())

#args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Rscript --vanilla ~/bin/loh.R normal.het.baf.txt tumor.ave.baf.txt outName", call.=FALSE)
}
library(ExomeCNV)
# object construction
normal = read.delim(args[1], header=TRUE)
tumor = read.delim(args[2], header=TRUE)
eLOH = LOH.analyze(normal, tumor, alpha=0.05, method="two.sample.fisher")
eLOH <- subset(x = eLOH, tumor.coverage > 4)
mLOH=multi.LOH.analyze(normal, tumor, all.loh.ls=list(eLOH), min.spec = 0.9999, method="two.sample.fisher")
mLOH <- subset(x = mLOH, subset = LOH)
pdf(file = paste(args[3], ".loh.pdf", sep = ""))
do.plot.loh(mLOH, normal, tumor, "two.sample.fisher", plot.style="dev")
dev.off()  
chr <- mLOH$chr
start <-mLOH$position.start
end <-mLOH$position.end
score <- -10 * log10(mLOH$pval) 
bLOH <- data.frame(chr, start, end, name = ".",  score, strand = ".")
bLOH$strand <- ifelse(bLOH$start <= bLOH$end, "+", "-")
n <- bLOH$strand=="-"
x <- bLOH$start[n]
bLOH$start[n] <- bLOH$end[n]
bLOH$end[n] <- x
write.table(bLOH, file=paste(args[3], ".loh.bed", sep = ""), row.names = FALSE, col.names = FALSE, quote=FALSE, sep = "\t")

