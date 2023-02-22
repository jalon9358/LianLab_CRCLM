#### spatial feature plot 20220416, by Long Jie ####
myspatial_featureplot=function(data,filename,genename,pt.size.factor=1.6,image.alpha=1){
  marker_no=length(genename)
  featureplot_no=marker_no/6
  featureplot_count=1
  while(featureplot_no>1){
    pdf(paste0(filename,"_","spatialFeature_",as.character(featureplot_count),".pdf"),width = 20, height = 10.3)
    g=SpatialFeaturePlot(object = data,features = genename[(6*(featureplot_count-1)+1):(6*featureplot_count)],
                         pt.size.factor=pt.size.factor,ncol = 3,
                         image.alpha = image.alpha)
    print(g)
    dev.off()
    png(paste0(filename,"_","spatialFeature_",as.character(featureplot_count),".png"),width = 3000, height = 2400,type = "cairo",res = 300)
    g=SpatialFeaturePlot(object = data,features = genename[(6*(featureplot_count-1)+1):(6*featureplot_count)],
                         pt.size.factor=pt.size.factor,ncol=3,image.alpha = image.alpha)
    print(g)
    dev.off()
    featureplot_no=featureplot_no-1
    featureplot_count=featureplot_count+1
  }
  if(featureplot_no>0){
    pdf(paste0(filename,"_","spatialFeature_",as.character(featureplot_count),".pdf"),
        width = 6.7*ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,                                                                                                 (marker_no-(6*(featureplot_count-1)))), 
        height = 5.1*ceiling((marker_no-(6*(featureplot_count-1)))/3))
    g=SpatialFeaturePlot(object = data, features = genename[(6*(featureplot_count-1)+1):marker_no],pt.size.factor=pt.size.factor,
                  ncol=ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,
                              (marker_no-(6*(featureplot_count-1)))),
                  image.alpha = image.alpha)
    print(g)
    dev.off()
    png(paste0(filename,"_","spatialFeature_",as.character(featureplot_count),".png"),
        width = 1000*ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,
                            (marker_no-(6*(featureplot_count-1)))), 
        height = 800*ceiling((marker_no-(6*(featureplot_count-1)))/3),
        type = "cairo",res = 300)
    g=SpatialFeaturePlot(object = data, features = genename[(6*(featureplot_count-1)+1):marker_no],
                  pt.size.factor=pt.size.factor,
                  ncol=ifelse((marker_no-(6*(featureplot_count-1)))>=3,3,
                              (marker_no-(6*(featureplot_count-1)))),
                  image.alpha = image.alpha)
    print(g)
    dev.off()
  }
} 
