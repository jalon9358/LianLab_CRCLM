##please set the choose coloum as factor----
percentplot_choose=function(seurat_object,filename,width=6,height=4,
                            legend.key.size=4,sample_name=NULL,flip=F,angle=45,
                            colors=colors2,choose=NULL){
  num <- table(seurat_object@active.ident)
  if(is.null(sample_name)){
    sample_name=levels(seurat_object@meta.data[,choose])
  }
  #sample_name=c("Normal","PE")
  frequency_matrix=data.frame(names(num))
  cluster=levels(seurat_object)
  for (i in sample_name) {
    cell_name=row.names(seurat_object@meta.data[seurat_object@meta.data[,choose]==i,])
    subset_data=subset(seurat_object, cells = cell_name)
    frequency=table(subset_data@active.ident)/ncol(subset_data)
    frequency=as.data.frame(frequency)
    zero=setdiff(levels(seurat_object@active.ident),
                 names(table(subset_data@active.ident)))
    zero=as.data.frame(zero)
    zero$Freq=c(rep(0,nrow(zero)))
    names(zero)[1]="Var1"
    frequency=rbind(frequency,zero)
    row.names(frequency)=frequency$Var1
    frequency=frequency[cluster,]
    frequency_matrix[,i]=frequency$Freq
  }
  names(frequency_matrix)=c("cluster_name",sample_name)
  row.names(frequency_matrix)=frequency_matrix[,1]
  write.csv(frequency_matrix,paste0(filename,"_frequency_matrix.csv"),row.names = F)
  library(reshape2)
  data_rownames <- row.names(frequency_matrix)
  data_colnames <- colnames(frequency_matrix)
  frequency_matrix$cluster <- data_rownames
  frequency_matrix=frequency_matrix[,-1]
  frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
  frequency_matrix_m <- melt(frequency_matrix, id.vars=c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
  colnames(frequency_matrix_m)=c("Cluster","Sample","Frequency","count","Percent")
  #colorm=cbind(names(num),colors[names(num)])
  colorm=cbind(names(num),colors[1:length(names(num))])
  row.names(colorm)=colorm[,1]
  color_percent=as.character(colorm[,2])
  frequency_matrix_m2=as.data.frame(frequency_matrix_m)
  p <- ggplot(frequency_matrix_m2, aes(x=Sample, y=Percent, group=Cluster)) +
    geom_bar(stat="identity", color="white",position="fill", aes(fill=Cluster))+
    scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = angle),
          legend.key.size=unit(legend.key.size,'mm'))
  if(flip==TRUE){
    p=p+coord_flip()
  }
  pdf(paste0("percent_plot_",filename,".pdf"),width =width, height = height*2)
  print(p)
  dev.off()
}