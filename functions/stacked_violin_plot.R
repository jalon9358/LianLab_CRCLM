#****Long Jie,2020-12-02****#
#This function is used to plot the stacked violinplot of single cell RNAseq data 
#or any seurat object for specified genes expressed on celltype clusters.

#two main argument must supply
#@gene:a vector of the genes name to plot
#@seurat_object:the seurat object of your data

#Optinal arguments
#@text.size:numeric,the text size of gene names
#@flip:logical,set the direction of coordinate axis
#@limits.max:numeric,the max legend value to regulate the color levels

#example shows below
###@@@@@@@####
#setwd("/Users/jalon/Desktop/")#set the export work directory
#gene=c("Cd44","Sell","Foxp3","Cd160","Klrg1")
#stacked_violin_plot(gene,seurat_object = T.subset2,text.size = 10,flip = F,
#                 filename = "myviolinplot",width = 12,height = 10,limits.max = 9,Mean=T,col=colors)
###@@@@@@@###
#library packages----
library(Seurat)
library(Rmisc)
library(ggplot2)

####**** function part ****####
stacked_violin_plot=function(gene,seurat_object,cluster=NULL,limits.max=7,
         width=13,height=10.3,flip=T,filename="",text.size=10,Mean=T,
         col=colors){
  if(length(cluster)>0){
    seurat_object=subset(seurat_object,idents = cluster)
  }
  ave_expression=AverageExpression(seurat_object,assays = "RNA")$RNA
  ave_expression=log2(ave_expression+1)
  data_matrix=seurat_object@assays$RNA@data
  plot.list=list()
  #g=gene[5]
  no=1
  gene=c(gene[1],gene,gene[length(gene)])
  for (g in gene) {
      ave_gene_choose=ave_expression[which(rownames(ave_expression)==g),]
      data_matrix_choose=data_matrix[which(rownames(data_matrix)==g),]
      #data_matrix_choose=as.data.frame(data_matrix_choose)
      df=data.frame(expression=as.numeric(data_matrix_choose),cluster=as.character(seurat_object@active.ident))
      mean=vector()
      ave_gene_choose=as.data.frame(ave_gene_choose)
      if(Mean==T){
        for (i in df$cluster) {
          mean=c(mean,ave_gene_choose[i,])
        }
        df[,"mean"]=as.data.frame(mean)
        df=as.data.frame(df)
      }
    
  if(!flip){
    df[,"cluster"]=factor(df$cluster,levels=levels(seurat_object))
    if(Mean==F){
      p <- ggplot(df, aes(x=cluster, y=expression, fill= cluster, color=cluster))+
        geom_violin(scale="width") +
        labs(title=paste(g), y ="Expression", x="Cluster")+
        #theme_classic() +
        scale_fill_manual(values = col)+
        scale_color_manual(values = col)+
        theme(axis.title.y =  element_blank())+
        #theme(axis.ticks.y =  element_blank())+
        #theme(axis.line.y =   element_blank())+
        #theme(axis.text.y =   element_blank())+
        theme(axis.title.x = element_blank())+
        theme(legend.position="none" )
    }else{
      p <- ggplot(df, aes(x=cluster, y=expression, fill= mean, color=mean))+
        geom_violin(scale="width") +
        labs(title=paste(g), y ="Expression", x="Cluster")+
        #theme_classic() +
        scale_color_gradientn(colors = c("#FFFF00", "#FFD000","#FF0000","#360101"),
                              limits=c(0,limits.max))+
        scale_fill_gradientn(colors = c("#FFFF00", "#FFD000","#FF0000","#360101"),
                             limits=c(0,limits.max))+
        theme(axis.title.y =  element_blank())+
        #theme(axis.ticks.y =  element_blank())+
        #theme(axis.line.y =   element_blank())+
        #theme(axis.text.y =   element_blank())+
        theme(axis.title.x = element_blank())+
        theme(legend.position="none" )
    }
    
    if(no!=length(gene)){
      p<-p+
        theme( axis.line.x=element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
    }else{
      p<-p+
        theme(axis.text.x = element_text(size = 10,vjust = 0.5,face = "bold",color = "black"))
    }
    #p<-p+theme(plot.title = element_text(size=text.size,face="bold",hjust = 0.5))
    p=p+theme(panel.border = element_rect(fill = "NA",size = 0.5,color = "black"),
              panel.background = element_blank(),
              panel.grid = element_blank())
    if(no==1){
      plot.margin=unit(c(0.2, 0.5, 0.2, 0.5), "cm")
      p=p+theme(legend.title = element_text(size = 9,face = "bold"), 
            legend.text  = element_text(size = 9,face = "bold"),
            legend.key.size = unit(1, "lines"))
      legend.position="none"
    }else if(no==length(gene)){
      plot.margin=unit(c(-0.3, 0.5, 0.2, 0.5), "cm")
      legend.position="none"
    }else{
      plot.margin=unit(c(-0.74, 0.5, 0, 0.5), "cm")
      legend.position="none"
    }
    p=p+xlab("") + ylab(g) + ggtitle("") +
      theme(legend.position = legend.position,
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            #axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_text(size=text.size,face="bold",hjust = 0.5),
            plot.margin = plot.margin )
    
    if(length(plot.list)==0){
      plot.list=list(p)
    }else{
      plot.list=c(plot.list,list(p))
    }
    no=no+1
  }else{
    ####******########
    df$cluster=factor(df$cluster,levels=rev(levels(seurat_object)))
    if(Mean==F){
      p <- ggplot(df, aes(x=cluster, y=expression, fill= cluster, color=cluster))+
        geom_violin(scale="width") +
        labs(title=paste(g), y ="Expression", x="Cluster")+
        #theme_classic() +
        scale_fill_manual(values = col)+
        scale_color_manual(values = col)+
        theme(axis.title.y =  element_blank())+
        #theme(axis.ticks.y =  element_blank())+
        #theme(axis.line.y =   element_blank())+
        #theme(axis.text.y =   element_blank())+
        theme(axis.title.x = element_blank())+
        theme(legend.position="right")
    }else{
      p <- ggplot(df, aes(x=cluster, y=expression, fill= mean, color=mean))+
        geom_violin(scale="width") +
        labs(title=paste(g), y ="Expression", x="Cluster")+
        #theme_classic() +
        scale_color_gradientn(colors = c("#FFFF00", "#FFD000","#FF0000","#360101"),
                              limits=c(0,limits.max))+
        scale_fill_gradientn(colors = c("#FFFF00", "#FFD000","#FF0000","#360101"),
                             limits=c(0,limits.max))+
        theme(axis.title.y =  element_blank())+
        #theme(axis.ticks.y =  element_blank())+
        #theme(axis.line.y =   element_blank())+
        #theme(axis.text.y =   element_blank())+
        theme(axis.title.x = element_blank())+
        theme(legend.position="right")
    }

    
    if(no!=length(gene)){
      p<-p+
        theme( axis.line.x=element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
    }else{
      p<-p+
        theme(axis.text.x = element_text(size = 10,vjust = 0.2,face = "bold",color = "black"))
    }
    p<-p+theme(plot.title = element_text(size=4,face="bold",hjust = 0.5,color = "black"))
    p=p+theme(panel.border = element_rect(fill = "NA",size = 0.5,color = "black"),
              panel.background = element_blank(),
              panel.grid = element_blank())
    if(no==1){
      plot.margin=unit(c(0.5, 0.5, 0.5, 0.2), "cm")
      p=p+theme(legend.title = element_text(size = 9,face = "bold"), 
                legend.text  = element_text(size = 9,face = "bold"),
                legend.key.size = unit(0.5, "lines"))
      #legend.position="left"
    }else if(no==length(gene)){
      plot.margin=unit(c(0.5, 0.2, 0.5, 0.3), "cm")
      #legend.position="none"
    }else{
      plot.margin=unit(c(0.5, 0, 0.5, -0.11), "cm")
      #legend.position="none"
    }
    if(no==1){
      p=p+xlab("") + ylab("") +
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              #axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              #axis.ticks.y = element_blank(),
              #axis.title.x = element_text(size=10,face="bold",hjust = 0.5),
              plot.title = element_text(colour = "black", face = "bold", 
                                        size = text.size, vjust = 0.2),
              axis.text.y = element_text(size=10,face="bold",hjust = 1,color = "black"),
              plot.margin = plot.margin )+coord_flip()
    }else if(no==length(gene)){
      p=p+xlab("") + ylab("") +
        theme(legend.position = "right",
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              #axis.title.x = element_text(size=10,face="bold",hjust = 0.5),
              plot.title = element_text(colour = "black", face = "bold", 
                                        size = text.size, vjust = 0.2),
              #axis.text.y = element_text(size=10,face="bold",hjust = 1,color = "black"),
              plot.margin = plot.margin )+coord_flip()
    }else{
      p=p+xlab("") + ylab("") +
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              #axis.title.x = element_text(size=10,face="bold",hjust = 0.5),
              plot.title = element_text(colour = "black", face = "bold", 
                                        size = text.size, vjust = 0.2),
              #axis.text.y = element_text(size=10,face="bold",hjust = 1,color = "black"),
              plot.margin = plot.margin )+coord_flip()
    }
    
    
    if(length(plot.list)==0){
      plot.list=list(p)
    }else{
      plot.list=c(plot.list,list(p))
    }
    no=no+1
  }
  }
  pdf(paste0(filename,"_vln_manual",".pdf"),width = width,height = height)
  if(flip){
    m=multiplot(plotlist = plot.list,cols = length(gene))
  }else{
    m=multiplot(plotlist = plot.list,cols=1)
  }
  dev.off()
}

