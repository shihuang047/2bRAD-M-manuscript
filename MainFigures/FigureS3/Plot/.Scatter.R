#--------------------------------------------------------------------------------------------
##  Scatter-plot For P3 Un-HomeCar
#Underarm
Un <- read.table("Underarm_long_data.txt",header = T,sep="\t")
ggplot(data=Un,aes(x=(log(X2B)),y=(log(X16S)),color=Color))+
  geom_point(aes(size=0.1),alpha=0.5) +
  theme_bw()+
  geom_smooth(method="lm",level=0.6) +
  scale_x_continuous(limits = c(-10,0))+scale_y_continuous(limits = c(-10,0)) +
  theme(legend.position="none") + 
  facet_wrap(~Group,scales = "free")
#Underarm_Eu
Un_E <- read.table("Underarm_long_data_Eu.txt",header = T,sep="\t")
ggplot(data=Un_E,aes(x=(log(X2B)),y=(log(X16S)),color=Color))+
  geom_point(size=3,alpha=0.5) +
  theme_bw()+
  geom_smooth(method="lm",level=0.97) +
  scale_x_continuous(limits = c(-10,0))+scale_y_continuous(limits = c(-10,0)) +
  theme(legend.position="none") + 
  facet_wrap(~Group,scales = "free")
#Home_car_Eu_Ar
Hc_EA <- read.table("Home_car_long_data_Eu_Ar.txt",header = T,sep="\t")
ggplot(data=Hc_EA,aes(x=(log(X2B)),y=(log(X16S)),color=Color))+
  geom_point(size=3,alpha=0.5) +
  theme_bw()+
  geom_smooth(method="lm",level=0.97) +
  scale_x_continuous(limits = c(-10,0))+scale_y_continuous(limits = c(-10,0)) +
  theme(legend.position="none") + 
  facet_wrap(~Group,scales = "free")
#Forearm_16S_Eu
F_16S_E <- read.table("Forearm_long_data_16S_Eu.txt",header = T,sep="\t")
ggplot(data=F_16S_E,aes(x=(log(X2B)),y=(log(X16S)),color=Color))+
  geom_point(size=3,alpha=0.5) +
  theme_bw()+
  geom_smooth(method="lm",level=0.97) +
  scale_x_continuous(limits = c(-10,0))+scale_y_continuous(limits = c(-10,0)) +
  theme(legend.position="none") + 
  facet_wrap(~Group,scales = "free")
#Forearm_Shotgun_Eu
F_S_E <- read.table("Forearm_long_data_Shotgun_Eu.txt",header = T,sep="\t")
ggplot(data=F_S_E,aes(x=(log(X2B)),y=(log(mOTUs2)),color=Color))+
  geom_point(size=3,alpha=0.5) +
  theme_bw()+
  geom_smooth(method="lm",level=0.97) +
  scale_x_continuous(limits = c(-10,0))+scale_y_continuous(limits = c(-10,0)) +
  theme(legend.position="none") + 
  facet_wrap(~Group,scales = "free")