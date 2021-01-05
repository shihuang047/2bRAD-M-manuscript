#--------------------------------------------------------------------------------------------
##  Scatter-plot For P3 Deep sequence

####  Deep_shotgun-2B  ####
##  Reshape data
#Shotgun
a_Ds <- read.table("Deep_seq_shotgun.txt", header = T, sep = "\t")
a_Ds <- melt(a_Ds)
write.table (a_Ds, file = "C:\\Users\\luck\\Desktop\\P_G\\1 output_long_data\\output_Deep_seq_shotgun_long_data.txt", sep ="\t", row.names =F, col.names =T, quote =F)
#2B-species
a_D2B_s <- read.table("Deep_seq_2B_species.txt", header = T, sep = "\t")
a_D2B_s <- melt(a_D2B_s)
write.table (a_D2B_s, file = "C:\\Users\\luck\\Desktop\\P_G\\1 output_long_data\\output_Deep_seq_2B_species_long_data.txt", sep ="\t", row.names =F, col.names =T, quote =F)
##  Plot
Ds <- read.table("Deep_seq_long_data_shotgun.txt",header = T,sep="\t")
ggplot(data=Ds,aes(x=(log(X2B_s)),y=(log(shotgun)),color=Color))+
  geom_point(size=3.1,alpha=0.5) +
  theme_bw()+geom_smooth(method="lm",level=0.97) +
  scale_x_continuous(limits = c(-10,0))+scale_y_continuous(limits = c(-10,0)) +
  theme(legend.position="none") + 
  facet_wrap(~Group,scales = "free")

####  Deep_16S-2B  ####
##  Reshape data
#16S
a_D16s <- read.table("Deep_seq_16S.txt", header = T, sep = "\t")
a_D16s <- melt(a_D16s)
write.table (a_D16s, file = "C:\\Users\\luck\\Desktop\\P_G\\1 output_long_data\\output_Deep_seq_16_long_data.txt", sep ="\t", row.names =F, col.names =T, quote =F)
#2B
a_D2B_g <- read.table("Deep_seq_2B_genus.txt", header = T, sep = "\t")
a_D2B_g <- melt(a_D2B_g)
write.table (a_D2B_g, file = "C:\\Users\\luck\\Desktop\\P_G\\1 output_long_data\\output_Deep_seq_2B_genus_long_data.txt", sep ="\t", row.names =F, col.names =T, quote =F)
##  Plot
D16S <- read.table("Deep_seq_long_data_16S.txt",header = T,sep="\t")
ggplot(data=D16S,aes(x=(log(X2B_g)),y=(log(X16S)),color=Color))+geom_point(size=3.1,alpha=0.5) +
  theme_bw()+geom_smooth(method="lm",level=0.97) +
  scale_x_continuous(limits = c(-10,0))+scale_y_continuous(limits = c(-10,0)) +
  theme(legend.position="none") + 
  facet_wrap(~Group,scales = "free")