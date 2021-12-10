## Fig3
library(ggplot2)
msa <- read.table("Plot_MSA_Distribution.txt", header = T, sep = "\t",row.names = 1)
#Plot
ggplot(msa, aes(matrix, abd,fill= species))+
  geom_bar(stat='identity',position='fill',width=0.8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  scale_fill_manual(values = colours)