## FigS9
library(pheatmap)
data<-read.table("Plot_Un_HC.txt",header=TRUE,sep="\t",row.names=1)
pheatmap(t(data),kmeans_k=NA,cellwidth=12,cellheight=12,
         fontsize_row=12,fontsize_col=12,
         filename="Fungi_-64.pdf",cluster_rows = FALSE,color = colorRampPalette(c( "white", "darkgoldenrod1", "red"))(20),
         angle_col = "315")