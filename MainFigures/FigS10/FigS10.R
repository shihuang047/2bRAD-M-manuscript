## FigS10
library(UpSetR)
library(RColorBrewer)

venn_dat <- read.table(args[1], header = T, sep = '\t', stringsAsFactors = F, check.names = F)

mycolors<- brewer.pal(5,"Set1")

pdf(paste(args[2],".Venn.pdf",sep=""), width=14, height=7, onefile=F)
upset(venn_dat, nsets = 9, number.angles = 30, point.size = 2, line.size = 1, sets.bar.color=mycolors, order.by = c("freq", "degree"))
dev.off()