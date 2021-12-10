## FigS6
library(ggplot2)
aa <- read.table("plot_rarefaction.txt",header = T,sep="\t")
ggplot(data=aa, aes(x=xa,y=ya,colour= SampleID))+
  geom_point(size = 3, alpha = 0.5)+theme_bw()+
  geom_smooth(level = 0.9, method="loess", formula = y ~ log(x))+
  facet_wrap(spp~type,scales = "free")