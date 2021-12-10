## FigS10
#!/usr/bin/env Rscript
args=commandArgs(T);
if (length(args) != 2 ){
	cat("\033[0;32m
	USAGE:
		Rscript UpsetR-venn.v1.R data Figuer
	AUTHOR:  Zhang Rongchao
	VERSION: 2019.03.16
	\033[0m\n"
	);
	q();
}

library(UpSetR)
library(RColorBrewer)

venn_dat <- read.table(args[1], header = T, sep = '\t', stringsAsFactors = F, check.names = F)

mycolors<- brewer.pal(5,"Set1")

pdf(paste(args[2],".Venn.pdf",sep=""), width=14, height=7, onefile=F)
upset(venn_dat, nsets = 9, number.angles = 30, point.size = 2, line.size = 1, sets.bar.color=mycolors, order.by = c("freq", "degree"))
dev.off()

