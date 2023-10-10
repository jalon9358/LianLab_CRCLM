#****Long Jie,2021-03-03****#
#This function is used to add the score of own genelist to the seurat object
#Two main argument must supply
#@genelist:the genes used to caculate the score
#@seurat_object:it also support a seurat_object instead of expression matrix
#example as belows:
##mymodule_score(seurat_object = myobject,genelist=mygenelist,cutoff=0.2,filename="mymodule")

library(tidyverse)
library(ggforce)
mymodule_score=function(seurat_object,genelist,cutoff=0.1,width.ppi=6,height.ppi=5,cutline=T,
                        filename="",pt.size=0.5){
  choose=subset(seurat_object,features=genelist)
  df=as.data.frame(choose@assays$RNA@scale.data)
  score=apply(df, 2, mean)
  score=log1p(score)
  
  pdf(paste0(filename,"_density_score.pdf"),width = width.ppi*1,height = height.ppi)
  base::plot(density(score))
  abline(v=cutoff)
  dev.off()
  
  data=GetAssayData(object = seurat_object, slot = "data")
  score=as.data.frame(score)
  score.df=score
  score=t(score)
  data=rbind(data,score)
  choose=seurat_object
  choose@assays$RNA@data=as.sparse(data)
  
  p2=ggplot(score.df,aes(x="",y=score))+
    geom_sina(color = "hotpink2", size = 0.01)+
    geom_boxplot(width = 0.1,color="hotpink2",
                 position = position_dodge(0.1))+theme_bw()+labs(y="Score")+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
    NoLegend()
  
  if(cutline){
    p2=p2+geom_hline(yintercept = cutoff,col="royalblue1",size=0.5,lty=2)
  }
  
  positive.per=length(which(score>cutoff))/length(score)
  p2=p2+labs(title = paste0("positive percentage is ",100*round(positive.per,4),"%"))+
    theme(plot.title = element_text(hjust = 0.5))
  
  pdf(paste0(filename,"_score_violin.pdf"),width = width.ppi*1,height = height.ppi)
  print(p2)
  dev.off()
      
  choose$score=t(score)
  choose$score_cluster=ifelse(choose$score>cutoff,"Positive","Negative")
  choose$score_cluster=factor(choose$score_cluster,levels = c("Positive","Negative"))
      
  pdf(paste0(filename,"_score_cluster_tsne.pdf"),width = width.ppi*1.3,height = height.ppi)
  p3=DimPlot(object = choose, reduction = 'tsne',label = F,cols = colors_crc,group.by = "score_cluster",
                     pt.size = pt.size)
  print(p3)
  dev.off()

  choose@assays$RNA@data=choose@assays$RNA@data[1:(nrow(choose)-1),]
  return(choose)
}
  