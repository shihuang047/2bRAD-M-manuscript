#--------------------------------------------------------------------------------------------
##    Heatmap For P3 Cancer
library(pheatmap)
library(RColorBrewer)
#Read table
data1<-read.table("Fungi_-6.txt",header=TRUE,sep="\t",row.names=1)
#Plot
pheatmap(t(data1),kmeans_k=NA,cellwidth=12,cellheight=12,
         fontsize_row=12,fontsize_col=12,
         filename="Fungi_-64.pdf",cluster_rows = FALSE,color = colorRampPalette(c( "white", "darkgoldenrod1", "red"))(20),
         angle_col = "315")#rainbow angle_col