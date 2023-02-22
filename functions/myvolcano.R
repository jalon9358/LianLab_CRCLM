library(ggrepel)
library(dplyr)
myvolcano=function(df,gene.plot=10,logfc.cutoff=0,p.cutoff=0.01,filename="",
                   width=6.5,height=5,text.size=3,pt.size=1,from.seurat=T,text.repel=T,
                   choose_genes=NULL){
  df$color = '#BEBEBE'
  df$group = 'Other'
  df$label = rownames(df)
  df$gene = rownames(df)
  if(from.seurat==F){
    colnames(df)[grep(pattern = "logFC",colnames(df))]="avg_log2FC"
    colnames(df)[grep(pattern = "P.Value",colnames(df))]="p_val_adj"
  }
  up = df %>% filter(., avg_log2FC > logfc.cutoff & (p_val_adj < p.cutoff))
  down = df %>% filter(.,avg_log2FC < -logfc.cutoff & (p_val_adj < p.cutoff))
  up=up[order(up$avg_log2FC,decreasing = T),]
  down=down[order(down$avg_log2FC,decreasing = F),]
  df[up$gene,'color'] = '#FF0000'
  df[up$gene,'group'] = 'up'
  df[down$gene,'color'] = '#0000FF'
  df[down$gene,'group'] = 'down'
  df$p_val_adj = -log10(df$p_val_adj)
  df = df[order(df$group,df$avg_log2FC),]
  DGE=df %>% group_by(group)%>%filter(.,group!="Other")
  #top = df %>% group_by(group) %>% top_n(n = gene.plot, wt = avg_log2FC)
  #top = top %>% filter(.,group !='Other')
  df[!df$gene %in% (c(up$gene[1:gene.plot],down$gene[1:gene.plot])),'label'] = ''
  if(!is.null(choose_genes)){
    df$label=ifelse(df$gene%in%choose_genes,df$gene,df$label)
  }
  
  pdf(file = paste(filename,'volcano_',gene.plot,'.pdf'),width = width,height = height)
  p <- ggplot(df, aes(avg_log2FC, p_val_adj,color = group)) + geom_point(size=pt.size)+
    labs(x='log2FC',y='-Log10 P')
  #labs(x='logFC',y='adjust p-value')
  p <- p + theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text = element_text(size = 16),axis.title = element_text(size = 16),
          axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
          axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
          axis.line = element_line(colour = 'black',size = 0.8),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = 'grey',size=0.6,linetype = 2),
          #panel.grid.major =element_blank(),
          axis.ticks = element_line(colour = 'black',size = 0.8),
          legend.title = element_blank(),legend.text = element_text(size = 14),
          legend.position = 'top',legend.key.size = unit(1.5, 'lines'))
  p = p + scale_color_manual(values = setNames(df$color, df$group))
  
  #p = p + geom_label_repel(aes(label = label,colour = 'red'),box.padding = 0.35,point.padding = 0.5,segment.color = 'grey50')
  if(text.repel){
    p = p + geom_text_repel(aes(label = label),color= 'black',size = text.size,segment.color = 'grey70')
  }
  #p = p + geom_hline(yintercept=0, size = 0.8, color = "grey")
  #p = p + geom_vline(xintercept=0, size = 0.8, color = "grey")
  #p = p + geom_abline(intercept = 0, slope = 1, color="grey", 
  #                    size=0.8)
  #p = p + guides(colour = guide_legend(override.aes = list(size = 5),nrow=3))
  p = p + theme(legend.key=element_rect(fill=NA))
  print(p)
  dev.off()
  up_gene=DGE[DGE$avg_log2FC>0,]$gene
  down_gene=DGE[DGE$avg_log2FC<0,]$gene
  export_gene=list(up=up_gene,down=down_gene)
  return(export_gene)
}
