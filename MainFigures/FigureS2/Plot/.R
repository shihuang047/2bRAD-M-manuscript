aa <- read.table("rarefraction_ZM_0717.txt",header = T,sep="\t")

ggplot(data=aa, aes(x=xa,y=ya,colour= SampleID))+
    geom_point(size = 1.5, alpha = 0.5)+theme_bw()+
    geom_smooth(level = 0.9, method="loess", formula = y ~ log(x))+
    facet_grid(spp~type,scales = "free")#+  scale_y_log10()