## Fig3
library(ggplot2)
cas <- read.table("plot_MOCK_CAS.txt", header = T, sep = "\t",row.names = 1)
#Plot
ggplot(cas, aes(Group, Abundance,fill= Species))+geom_bar(stat='identity',position='fill',width = 0.7)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  scale_fill_manual(values = colours)+
  geom_text(aes(label=Dist, y=Abundance+0),angle = 90)