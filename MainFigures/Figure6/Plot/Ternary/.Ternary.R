aa <- read.table("suboptimal_Ternary_cancer_RF_HS_35.txt",header = T, sep="\t")
ggtern(aa,aes(CA, CIN, Nor)) +
  theme_rgbw() +
  geom_mask() +
  geom_point(size=3,aes(shape=Group,fill=Group)) +
  scale_shape_manual(values=c(21,24,18)) +
  labs(title="Demonstration of Raster Annotation")