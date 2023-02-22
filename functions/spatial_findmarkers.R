#### myfindmarkers for spatial rna-seq data, 20220415 by Long Jie ####
spatial_findmarkers=function(seurat_object,filename,thresh.use=0.25,species="human",
                             all.markers=NULL,colors=colors2,
                       fill_colors=colorRampPalette(c("lightskyblue", "white", "red"))(100)){
  if(is.null(all.markers)){
    all.markers <- FindAllMarkers(object = seurat_object, only.pos = TRUE, min.pct = 0.25, 
                                  thresh.use = thresh.use)
  }
  
  common_markers=c()
  if(length(grep(pattern = "TotalSeqC",rownames(seurat_object)))>0){
    antibody_name=rownames(seurat_object)[grep(pattern = "TotalSeqC",rownames(seurat_object))]
    common_markers=c(common_markers,antibody_name)
  }
  non_common=setdiff(unique(as.character(all.markers$gene)),common_markers)
  all.markers=subset(all.markers,gene%in%non_common)
  top5 <- all.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
  top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
  top50 <- all.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
  write.csv(top20,paste0(filename,"_clustertop20_gene.csv"),row.names = F)
  write.csv(top50,paste0(filename,"_clustertop50_gene.csv"),row.names = F)
  heatgene=unique(as.character(top20$gene))

  DefaultAssay(seurat_object)="SCT"
  seurat_object=ScaleData(seurat_object,features = rownames(seurat_object))
  
  ave_cluster=AverageExpression(seurat_object,assays = "SCT",slot = "data",return.seurat = T)
  pdf(paste0(filename,"_heatmap","_top20.pdf"),width =15, height = 40)
  f1=DoHeatmap(ave_cluster, features = heatgene,
               label = TRUE,size=3,assay = "SCT",lines.width = 1,
               group.colors = colors,draw.lines = F)+
    scale_fill_gradientn(colors = fill_colors)
  print(f1)
  dev.off()
  
  pdf(paste0(filename,"_heatmap","_width_top20.pdf"),width =16, height = 12)
  print(f1)
  dev.off()
  
  pdf(paste0(filename,"_heatmap","_top5.pdf"),width =15, height = 12)
  f2=DoHeatmap(object = ave_cluster, features = unique(as.character(top5$gene)),
               label = TRUE,size=3,group.colors = colors,draw.lines = F)+
    scale_fill_gradientn(colors = fill_colors)
  print(f2)
  dev.off()
  
  pdf(paste0(filename,"_heatmap","_width_top5.pdf"),width =16, height = 12)
  print(f2)
  dev.off()
  
  png(paste0(filename,"_heatmap","_top20.png"),width = 1500, height = 4000)
  print(f1)
  dev.off()
  
  png(paste0(filename,"_heatmap","_top5.png"),width = 750/1.5, height = 600/1.5)
  print(f2)
  dev.off()
  
  pdf(paste0(filename,"_heatmap","_nolegend.pdf"),width =20, height = 5.18)
  f3=DotPlot(object = seurat_object, features = rev(x = unique(as.character(top5$gene))), cols = c("blue","firebrick2"), 
             dot.scale = 5) + RotatedAxis()+NoLegend()
  print(f3)
  dev.off()
  
  pdf(paste0(filename,"_total_heatmap","_top5.pdf"),width =15, height = 12)
  f4=DoHeatmap(object = seurat_object, features = unique(as.character(top5$gene)),
               label = TRUE,size=3,group.colors = colors,draw.lines = T)+
    scale_fill_gradientn(colors = fill_colors)
  print(f4)
  dev.off()
  
  pdf(paste0(filename,"_total_heatmap","_top20.pdf"),width =15, height = 40)
  f5=DoHeatmap(object = seurat_object, features = unique(as.character(top20$gene)),
               label = TRUE,size=3,group.colors = colors,draw.lines = T)+
    scale_fill_gradientn(colors = fill_colors)
  print(f5)
  dev.off()
  
  export=list(top5=unique(as.character(top5$gene)),top20=unique(as.character(top20$gene)))
  return(export)
}

