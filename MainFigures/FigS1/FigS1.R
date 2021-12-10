## FigS1
library(ggplot2)
ph <- read.table("phylum.txt",header=T,sep="\t")
ggplot(ph,aes(x=enzyme,y=log10(all)))+geom_violin(aes(fill=enzyme))+
  labs(title="phylum")+geom_boxplot(width = 0.2)+theme_bw()+
  scale_y_continuous(breaks = c(0,2,5,10,15))