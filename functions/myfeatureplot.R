myfeatureplot=function(data,filename,genename,pt.size=0.5,min.cutoff="q5",
                       reduction="umap",max.cutoff="q99",
                       cols=c("lightblue1","lightgreen","darkorange","firebrick1"),
                       width=6.7,height=5.1){
  marker_no=length(genename)
  featureplot_no=marker_no/6
  featureplot_count=1
  while(featureplot_no>1){
    pdf(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".pdf"),width = 20, height = 10.3)
    g=FeaturePlot(object = data, reduction = reduction,features = genename[(6*(featureplot_count-1)+1):(6*featureplot_count)],cols = cols,pt.size = pt.size,ncol = 3,min.cutoff = min.cutoff,max.cutoff = max.cutoff,raster=F)
    print(g)
    dev.off()
    png(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".png"),width = 3000, height = 2400,type = "cairo",res = 300)
    g=FeaturePlot(object = data,reduction = reduction,features = genename[(6*(featureplot_count-1)+1):(6*featureplot_count)],cols = cols, pt.size = pt.size,ncol=3,min.cutoff = min.cutoff,max.cutoff = max.cutoff)
    print(g)
    dev.off()
    featureplot_no=featureplot_no-1
    featureplot_count=featureplot_count+1
  }
  if(featureplot_no>0){
    pdf(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".pdf"),width = width*ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,(marker_no-(6*(featureplot_count-1)))), height = height*ceiling((marker_no-(6*(featureplot_count-1)))/3))
    g=FeaturePlot(object = data,reduction = reduction, features = genename[(6*(featureplot_count-1)+1):marker_no],cols = cols, min.cutoff = min.cutoff,
                  max.cutoff = max.cutoff,pt.size = pt.size,ncol=ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,(marker_no-(6*(featureplot_count-1)))),raster=F)
    print(g)
    dev.off()
    png(paste0("Featureplot_",filename,"_",as.character(featureplot_count),".png"),width = 1000*ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,(marker_no-(6*(featureplot_count-1)))), height = 800*ceiling((marker_no-(6*(featureplot_count-1)))/3),type = "cairo",res = 300)
    g=FeaturePlot(object = data,reduction = reduction, features = genename[(6*(featureplot_count-1)+1):marker_no],cols = cols, 
                  max.cutoff = max.cutoff,min.cutoff = min.cutoff,pt.size = pt.size,ncol=ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,(marker_no-(6*(featureplot_count-1)))))
    print(g)
    dev.off()
  }
}