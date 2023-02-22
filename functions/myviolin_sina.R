## violin sina plot, 20220613, by Long Jie ----
myviolin_sina=function(data.object,colors,violin_gene,compare=F,comparisons=NULL,
                       module.score=FALSE,max.width=0.8,method = "wilcox.test"){
  if(module.score){
    data.matrix=data.object@meta.data[,violin_gene]
  }else{
    data.matrix=data.object@assays$RNA@data[violin_gene,]
  }
  data.matrix=as.data.frame(data.matrix)
  colnames(data.matrix)=violin_gene
  data.matrix$cluster=data.object@active.ident
  p=ggplot(data.matrix,aes(cluster,data.matrix[,violin_gene],fill=cluster,color=cluster))+theme_bw()+
    scale_color_manual(values = colors)+
    scale_fill_manual(values = colors)+
    #labs(x="GSE50760")+
    geom_violin(scale = "width",alpha=0,size=0.5)+
    geom_sina(size = 0.5,maxwidth=max.width,scale="width")+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.background = element_blank(),axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
    labs(y=violin_gene)
  if(compare){
    p=p+stat_compare_means(comparisons = list(comparisons),
                           label = "p.format",method = method,bracket.size = 0.5,
                           label.y = max(data.matrix[,violin_gene]),
                           size=5,paired = F)
  }
  print(p)
}
