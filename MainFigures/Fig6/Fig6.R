## Fig6a-b
PM_Marker_Test.R

## Fig6c
library(ggtern)
aa <- read.table("plot_Ternary.txt",header = T, sep="\t")
ggtern(aa,aes(CA, CIN, Nor)) +
  theme_rgbw() +
  geom_mask() +
  geom_point(size=3,aes(shape=Group,fill=Group)) +
  scale_shape_manual(values=c(21,24,18)) +
  labs(title="Demonstration of Raster Annotation")

## Fig6d
RF Model