## Fig4a-f
library(ggplot2)
Ds <- read.table("plot_scatter.txt",header = T,sep="\t")
ggplot(data=Ds,aes(x=(log10(X2B_g)),y=(log10(X16S)),color=Color))+
  geom_point(size=3.1,alpha=0.5) +
  theme_bw()+geom_smooth(method="lm",level=0.97) +
  scale_x_continuous(limits = c(-4.5,0))+scale_y_continuous(limits = c(-4.5,0)) +
  theme(legend.position="top") + 
  facet_wrap(~Group,scales = "free")

## Fig4g
library(ggplot2)
ggplot(data=aa, aes(x=xa,y=ya,colour= SampleID))+
  geom_point(size = 1.5, alpha = 0.5)+theme_bw()+
  geom_smooth(level = 0.9, method="loess", formula = y ~ log(x))+
  facet_grid(spp~type,scales = "free")+scale_x_continuous(breaks = c(0,5,10,13.5,15))
