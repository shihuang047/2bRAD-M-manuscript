## Fig2A
library(ggplot2)
ggplot(data=aa,aes(x=abs(Truth_tax_abd),y=abs(X2B),color=Color))+geom_point(aes(size=1),alpha=0.5) +
    theme_bw()+geom_smooth(method="lm",level=0.97) +
    scale_x_continuous(limits = c(0,0.04))+scale_y_continuous(limits = c(0,0.04)) +
    theme(legend.position="none")+facet_wrap(~Facet,scales = "free")