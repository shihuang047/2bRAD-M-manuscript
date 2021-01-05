#----------------------------------------------------------------------------------------
##    barplot--Distribution For P3 MOCK

#### CAS  ####
library(ggplot2)
library(RColorBrewer)
##  reshape data
a <- read.table("MOCK_CAS_distribution_noother.txt", header = T, sep = "\t")
a <- melt(a)
write.table (a, file = "C:\\Users\\luck\\Desktop\\long_data.txt", sep ="\t", row.names =F, col.names =T, quote =F)
##  Plot
#Read table
cas <- read.table("MOCK_CAS_distribution_noother_long_data.txt", header = T, sep = "\t",row.names = 1)
##  Change colour
#display.brewer.all(type = "all")
colours <- c(brewer.pal(7, "Set1"),brewer.pal(6, "Set2"),brewer.pal(8, "Accent")[c(5,7)],
             brewer.pal(8, "Dark2")[c(1,5)],brewer.pal(9, "YlOrRd")[c(3,6)],
             brewer.pal(12, "Paired")[c(2,12)])
#Plot
ggplot(cas, aes(Group, Abundance,fill= Species))+geom_bar(stat='identity',position='fill',width = 0.7)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  scale_fill_manual(values = colours)+
  geom_text(aes(label=Dist, y=Abundance+0),angle = 90)


#### MSA  ####
library(ggplot2)
library(RColorBrewer)
#Read table
msa <- read.table("MOCK_MSA_long_data_noother.txt", header = T, sep = "\t",row.names = 1)
##  Change order
msa$Species <- fct_inorder(msa$Species)#reverse order - distribution order by species name
##  Change colour
#display.brewer.all(type = "all")
colours <- c(brewer.pal(7, "Set1"),brewer.pal(6, "Set2"),brewer.pal(8, "Accent")[c(5,7)],
             brewer.pal(8, "Dark2")[c(1,5)],brewer.pal(9, "YlOrRd")[c(3,6)],
             brewer.pal(12, "Paired")[c(2,12)])
#Plot
ggplot(msa, aes(Group, Abundance,fill= Species))+
  geom_bar(stat='identity',position='fill',width=0.8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  scale_fill_manual(values = colours)+
  geom_text(aes(label=Dist, y=Abundance+0.75),angle = 90)