#example 
# setwd("/Users/jalon/Desktop/help_analysis/LL_AA/20210808_analysis/cpdb_merge_cluster/out")
# interaction_matrix=interaction_calculation(sig.matrix=sig.means)
# bubble_matrix=interaction_bubble(pvalue.path = "pvalues.txt",means.path = "means.txt",filename = "example",
#                                  pattern = "^M.*\\|T.*_CD8.*")
# 
# interaction_circlize(interaction_matrix$sig.no,pattern = "CD8",
#                      color = colors_aa,source.object = AA.combined.filter,filename = "example")
# 
# 
# mycluster_network(interaction_matrix$sig.no,pattern = "CD8",color = colors_aa,
#                   raw.object = AA.combined.filter,
#                   filename = "example_network",layout.factor = 1,
#                   edge.label = F,cutoff = 40)
# mypvals=order_LR_pairs(mypvals)
# mymeans=order_LR_pairs(mymeans)
# sig.means=order_LR_pairs(sig.means,sig.matrix = T)


##interaction no and summary score calculation function----
library(reshape2)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(igraph)
library(dplyr)
interaction_calculation=function(sig.matrix,means.cutoff=NULL){
  #sig.means=read.delim(sig.means.path,
   #                    check.names = F) #import significance matrix,default pvalue cutoff is 0.05
  sig.means=sig.matrix
  sig.means[is.na(sig.means)]=0 #set NA to 0
  sig.means=sig.means%>%dplyr::select("interacting_pair", #choose pair data
                                      grep(pattern = "\\|",colnames(sig.means),value = T))
  rownames(sig.means)=sig.means$interacting_pair
  sig.means=sig.means[,-1] #delete pair column
  if(!is.null(means.cutoff)){
    sig.means[sig.means<means.cutoff]=0
  }
  sig.summary=colSums(sig.means)
  sig.means.interno=colSums(sig.means>0)
  sig.no.matrix=data.frame(cluster1=str_split_fixed(names(sig.means.interno),"\\|",2)[,1],
                           cluster2=str_split_fixed(names(sig.means.interno),"\\|",2)[,2],
                           value=as.numeric(sig.means.interno))
  sig.summary.matrix=data.frame(cluster1=str_split_fixed(names(sig.summary),"\\|",2)[,1],
                                cluster2=str_split_fixed(names(sig.summary),"\\|",2)[,2],
                                value=as.numeric(sig.summary))
  sig.no.matrix.w<-reshape2::dcast(sig.no.matrix,cluster1~cluster2,value.var = 'value')
  rownames(sig.no.matrix.w)=sig.no.matrix.w$cluster1
  sig.no.matrix.w=sig.no.matrix.w[,-1]
  for(i in 1:nrow(sig.no.matrix.w)){
    for (j in 1:nrow(sig.no.matrix.w)) {
      if(i!=j){
        number=sig.no.matrix.w[i,j]+sig.no.matrix.w[j,i]
        sig.no.matrix.w[i,j]=number
        sig.no.matrix.w[j,i]=number
      }
    }
  }
  sig.summary.matrix.w<-reshape2::dcast(sig.summary.matrix,cluster1~cluster2,value.var = 'value')
  rownames(sig.summary.matrix.w)=sig.summary.matrix.w$cluster1
  sig.summary.matrix.w=sig.summary.matrix.w[,-1]
  for(i in 1:nrow(sig.summary.matrix.w)){
    for (j in 1:nrow(sig.summary.matrix.w)) {
      if(i!=j){
        number=sig.summary.matrix.w[i,j]+sig.summary.matrix.w[j,i]
        sig.summary.matrix.w[i,j]=number
        sig.summary.matrix.w[j,i]=number
      }
    }
  }
  list=list(sig.no=sig.no.matrix.w,sig.summary=sig.summary.matrix.w)
  return(list)
}


##cellphonedb bubble plot function part----
interaction_bubble=function(pvalue.path,means.path,pattern=NULL,
                            pvalue.cutoff=0.05,means.cutoff=0,
                            filename="",width=10,height=8,legend.midpoint=1,
                            size.factor=0.0001,scale.size=5,hjust=0.1,vjust=0.8){
  mypvals <- read.delim(pvalue.path, check.names = FALSE)
  mymeans <- read.delim(means.path, check.names = FALSE)
  mymeans %>%
    dplyr::select("interacting_pair",grep(pattern = pattern,
                                          colnames(mymeans),value = T))  %>%  
    reshape2::melt() -> meansdf
  colnames(meansdf)<- c("interacting_pair","CC","means")
  
  mypvals %>% 
    dplyr::select("interacting_pair",grep(pattern = pattern,
                                          colnames(mypvals),value = T))%>%  
    reshape2::melt()-> pvalsdf
  colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
  
  pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf <- merge(pvalsdf,meansdf,by = "joinlab")
  pdf(paste0("cellphone_bubble_",filename,".pdf"),width = width,height = height)
  p=pldf%>% filter(pvals<=pvalue.cutoff) %>%filter(means>=means.cutoff)%>%
    ggplot(aes(CC.x,interacting_pair.x) )+ 
    geom_point(aes(color=means,size=-log10(pvals+size.factor)) ) +
    scale_size_continuous(range = c(1,scale.size))+
    scale_color_gradient2(high="#AB1961",mid = "#FFD5B4",low ="#004E8F",midpoint = legend.midpoint)+ theme_bw()+ 
    theme(axis.text.x = element_text(angle = -90,hjust = hjust,vjust = vjust))
  print(p)
  dev.off()
  return(pldf)
}

##interaction circlize plot function----
interaction_circlize=function(matrix,pattern=NULL,color=colors,
                              source.object,filename="",width=12,
                              height=10,reduce=0.01,link.cutoff=0,
                              trackheight=c(0.03,0.05),transparency=0.2){
  if(!is.null(pattern)){
    matrix=subset(sig.matrix.w,subset=!grepl(pattern = pattern,rownames(matrix)),
                  select=grep(pattern = pattern,rownames(matrix)))
    matrix=t(matrix)
  }
  #colnames(sig.matrix.cir)=str_split_fixed(colnames(sig.matrix.cir),"_",2)[,1]
  grid.col=color
  names(grid.col)=levels(source.object)
  
  pdf(paste0("interaction_number_circlize_",filename,".pdf"),width = width,height = height)
  chordDiagram(as.matrix(matrix),reduce=reduce,link.visible = as.matrix(sig.matrix.cir)>=link.cutoff,
               annotationTrack = c("name", "grid"),annotationTrackHeight = trackheight,transparency = transparency,
               grid.col = grid.col)
  circos.clear()
  dev.off()
}

##igraph network cell interaction----
##network----
mycluster_network=function(matrix,pattern=NULL,cutoff=NULL,color=colors,raw.object,
                           max.vertex.size=20,max.edge.width=10,edge.label=TRUE,
                           layout.factor=1,vertex.label.cex=1,filename="defaultname",
                           width=10,height=8){
  matrix$cluster=rownames(matrix)
  matrix.long=melt(matrix)
  if(!is.null(pattern)){
    matrix.long=subset(matrix.long,subset=grepl(pattern = pattern,
                                                matrix.long$variable))
  }
  color1=color
  names(color1)=levels(raw.object)
  max=max(matrix.long$value)
  color2=colorRampPalette(brewer.pal(9, "Reds")[c(1,4,7)])(max) #将颜色分成多少份，取决于互作关系数目的最大值
  names(color2)=1:max #每一份颜色用对应的数字命名
  net <- graph_from_data_frame(matrix.long)
  matrix.width=reshape2::dcast(matrix.long,cluster~variable)
  rownames(matrix.width)=matrix.width$cluster
  matrix.width=matrix.width[,-1]
  ##set vertex size
  vsize=rowSums(matrix.width)
  vsize=vsize[names(V(net))]
  if(max(vsize)>max.vertex.size){
    vsize.factor=max(vsize)/max.vertex.size
    V(net)$size=vsize/vsize.factor
  }else{
    V(net)$size=vsize
  }
  #V(net)$size <- degree(net)*2+5  #节点大小与点中心度成正比，中心度即与该点相连的点的总数
  ##set edge width
  ewidth=E(net)$value
  if(max(ewidth)>max.edge.width){
    E(net)$width=ewidth*max.edge.width/max(ewidth)
  }else{
    E(net)$width=ewidth
  }
  E(net)$label.color <- "black" #连线标注的颜色
  V(net)$label.color <- 'black' #设置节点标记的颜色
  if(edge.label){
    E(net)$label <- E(net)$value #根据频次列设置边标签
  }else{
    E(net)$label=""
  }
  
  V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息
  E(net)$color <- color2[as.character(ifelse(E(net)$value > max,max,E(net)$value))] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
  #E(net)$arrow.size=0.3 #设置箭头大小
  #edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
  #调整节点位置的线条角度
  ##如果没有这两行代码，节点位置的圆圈是向右的
  #loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
  #igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  #l<-layout.fruchterman.reingold(net) #设置图的布局方式为弹簧式发散的布局
  #l <- layout_randomly(net)##修改布局
  #l=layout_on_sphere(net)
  #l=layout_with_fr(net)
  #l=layout_with_lgl(net)
  l=layout_with_drl(net)
  #l=layout_with_kk(net)
  #l=layout_with_graphopt(net)
  l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)#调整布局位置
  if(!is.null(cutoff)){
    net <- delete_edges(net, E(net)[value<cutoff])#删除边
  }
  #cut.off <- mean(sig.matrix.long.choose$value) #确定边的删除阈值
  pdf(paste0("interaction_network_",filename,".pdf"),width = width,height = height)
  plot(net,
       edge.arrow.size = 0, #连线不带箭头
       edge.curved = 0, #连线不弯曲
       vertex.frame.color = "black", #节点外框颜色
       layout = l*layout.factor,
       rescale=F,#使用设置位置
       vertex.label.cex = vertex.label.cex)#节点标注字体大
  dev.off()
}

###re-order ligand-receptor pairs-----
# order_LR_pairs=function(df,sig.matrix=FALSE){
#   out.df=data.frame()
#   if(sig.matrix){
#     metacol=12
#   }else{
#     metacol=11
#   }
#   for (i in 1:length(df$gene_a)) {
#     sub.df=df[i,]
#     if(sub.df$receptor_b=="False"){
#       if(sub.df$receptor_a=="True"){
#         oldnames=colnames(sub.df)
#         my_list=strsplit(oldnames[-c(1:metacol)],split = "\\|")
#         my_character=paste(sapply(my_list,'[[',2L),sapply(my_list,"[[",1L),sep = "|")
#         if(sig.matrix){
#           newnames=c(names(sub.df)[1:4],"gene_b","gene_a","secreted","receptor_b",
#                      "receptor_a","annotation_strategy","is_integrin","rank",my_character)
#         }else{
#           newnames=c(names(sub.df)[1:4],"gene_b","gene_a","secreted","receptor_b",
#                      "receptor_a","annotation_strategy","is_integrin",my_character)
#         }
#         sub.df=dplyr::select(sub.df,all_of(newnames))
#         names(sub.df)=oldnames
#         out.df=rbind(out.df,sub.df)
#       }
#     }else{
#       out.df=rbind(out.df,sub.df)
#     }
#   }
#   return(out.df)
# }

###re-order ligand-receptor pairs-----
order_LR_pairs=function(df,sig.matrix=FALSE){
  out.df=data.frame()
  if(sig.matrix){
    metacol=12
  }else{
    metacol=11
  }
  for (i in 1:length(df$gene_a)) {
    sub.df=df[i,]
    if(sub.df$receptor_b=="False"){
      if(sub.df$receptor_a=="True"){
        oldnames=colnames(sub.df)
        my_list=strsplit(oldnames[-c(1:metacol)],split = "\\|")
        my_character=paste(sapply(my_list,'[[',2L),sapply(my_list,"[[",1L),sep = "|")
        if(sig.matrix){
          newnames=c(names(sub.df)[c(1:2,4,3)],"gene_b","gene_a","secreted","receptor_b",
                     "receptor_a","annotation_strategy","is_integrin","rank",my_character)
        }else{
          newnames=c(names(sub.df)[c(1:2,4,3)],"gene_b","gene_a","secreted","receptor_b",
                     "receptor_a","annotation_strategy","is_integrin",my_character)
        }
        #sub.df=dplyr::select(sub.df,all_of(newnames))
        colnames(sub.df)=newnames
        #names(sub.df)[1:metacol]=oldnames[1:metacol]
        sub.df$interacting_pair=paste0(sub.df$gene_a,"_",sub.df$gene_b)
        sub.df=sub.df[,oldnames]
        out.df=rbind(out.df,sub.df)
      }else{
        out.df=rbind(out.df,sub.df)
      }
    }else{
      out.df=rbind(out.df,sub.df)
    }
  }
  return(out.df)
}

## 20220304 function ----
interaction_bubble_choose=function(pvalue.path,means.path,pattern=NULL,
                            pvalue.cutoff=0.05,means.cutoff=0,
                            filename="",width=10,height=8,legend.midpoint=1,
                            size.factor=0.0001,scale.size=5,hjust=0,vjust=0.5,choose.interaction=NULL){
  mypvals <- read.delim(pvalue.path, check.names = FALSE)
  mymeans <- read.delim(means.path, check.names = FALSE)
  mymeans %>%
    dplyr::select("interacting_pair",grep(pattern = pattern,
                                          colnames(mymeans),value = T))  %>%  
    reshape2::melt() -> meansdf
  colnames(meansdf)<- c("interacting_pair","CC","means")
  
  mypvals %>% 
    dplyr::select("interacting_pair",grep(pattern = pattern,
                                          colnames(mypvals),value = T))%>%  
    reshape2::melt()-> pvalsdf
  colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
  
  pvalsdf=subset(pvalsdf,subset=interacting_pair%in%choose.interaction)
  meansdf=subset(meansdf,subset=interacting_pair%in%choose.interaction)
  
  pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf <- merge(pvalsdf,meansdf,by = "joinlab")
  pdf(paste0("cellphone_bubble_",filename,".pdf"),width = width,height = height)
  p=pldf%>% filter(pvals<=pvalue.cutoff) %>%filter(means>=means.cutoff)%>%
    ggplot(aes(CC.x,interacting_pair.x) )+ 
    geom_point(aes(color=means,size=-log10(pvals+size.factor)) ) +
    scale_size_continuous(range = c(1,scale.size))+
    scale_color_gradient2(high="#AB1961",mid = "#FFD5B4",low ="#004E8F",midpoint = legend.midpoint)+ theme_bw()+ 
    theme(axis.text.x = element_text(angle = -90,hjust = hjust,vjust = vjust))
  print(p)
  dev.off()
  return(pldf)
}



