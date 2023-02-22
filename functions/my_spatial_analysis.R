## standard spatial analysis pretreatment pipeline,20220422 by Long Jie ##
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
source("E:/important/newest scripts/functions/source.R")
source("E:/important/newest scripts/functions/spatial_findmarkers.R")
source("E:/important/newest scripts/functions/myspatial_featureplot.R")
## function part ----
my_spatial_analysis=function(data.path,ndims=30,resolution=0.5,colors=colors_crc_spatial,
                             pt.size.factor = 1.4,width.ppi=6.5,
                             height.ppi=5,filename="defaultname"){
  cat("Importing datasets ... \n")
  spatial_object=Load10X_Spatial(data.path) ##load dataset
  spatial_object=subset(spatial_object,nFeature_Spatial>200)
  plot1 <- VlnPlot(spatial_object,features = "nCount_Spatial",pt.size = 0.1)+NoLegend()
  plot2 <- SpatialFeaturePlot(spatial_object,features = "nCount_Spatial")+
    theme(legend.position = "right")
  ggsave(paste0(filename,"_ncount_spatial.pdf"),wrap_plots(plot1, plot2),
         width = width.ppi,height = height.ppi)
  spatial_object=SCTransform(spatial_object, 
                             assay = "Spatial", verbose = FALSE) ##normalization
  pdf(paste0(filename,"_EPCAM_PTPRC.pdf"),width = width.ppi*1.8,height = height.ppi)
  p1=SpatialFeaturePlot(spatial_object, features = c("EPCAM", "PTPRC"),pt.size.factor = pt.size.factor)
  print(p1)
  dev.off()
  pdf(paste0(filename,"_ALB_PABPC1.pdf"),width = width.ppi*1.8,height = height.ppi)
  p2=SpatialFeaturePlot(spatial_object, features = c("ALB", "PABPC1"),pt.size.factor = pt.size.factor)
  print(p2)
  dev.off()
  ## dimension reduction ----
  spatial_object <- RunPCA(spatial_object, assay = "SCT", verbose = FALSE)
  pdf(paste0(filename,"_ElbowPlot.pdf"),width = width.ppi,height = height.ppi)
  p3=ElbowPlot(spatial_object,ndims = 50)
  print(p3)
  dev.off()
  spatial_object <- FindNeighbors(spatial_object, reduction = "pca", dims = 1:ndims)
  spatial_object <- FindClusters(spatial_object, verbose = FALSE,resolution = resolution)
  spatial_object <- RunUMAP(spatial_object, reduction = "pca", dims = 1:ndims)
  p4 <- DimPlot(spatial_object, reduction = "umap", label = F,cols = colors)
  p5 <- SpatialDimPlot(spatial_object,pt.size.factor = pt.size.factor,label = F, 
                       label.size = 3,cols = colors)
  pdf(paste0(filename,"_spatialdimplot.pdf"),width = width.ppi*1.8,height = height.ppi)
  p6=p4 + p5
  print(p6)
  dev.off()
  ## find marker genes ----
  cat("Caculating marker genes ... \n")
  spatial_markers=spatial_findmarkers(spatial_object,filename = filename,colors = colors)
  heatgene=unique(as.character(spatial_markers$top5))
  myspatial_featureplot(spatial_object,filename = filename,genename = heatgene,
                        pt.size.factor = pt.size.factor)
  
  crc_common_markers=c("COL1A1","S100A8","VGF",
                       "DES","LTB","DEFA6","SERF2",
                       "MUC4","FABP4","MCAM")
  functional.markers=c("SPP1","MMP2","CCL14","CXCL14","MKI67","CXCL8","CCL21",
                       "IFNG","IL10","TGFB1","IL1B","MCAM","NOTCH1","JAG1",
                       "NOTCH3","DLL4")
  immune.markers=c("SPP1","CXCL9","CD68","LYZ","MS4A1","CD79A","GNLY","GZMB",
                   "IL21","FOXP3","IGHM","LAMP3","CD14","CD3D","CD4","CD8A","IL2RA",
                   "JCHAIN","MZB1")
  myspatial_featureplot(spatial_object,filename = "functional_markers",genename = functional.markers,
                        pt.size.factor = pt.size.factor)
  myspatial_featureplot(spatial_object,filename = "immune_markers",genename = immune.markers,
                        pt.size.factor = pt.size.factor)
  myspatial_featureplot(spatial_object,filename = "crc_common_markers",genename = crc_common_markers,
                        pt.size.factor = pt.size.factor)
  #spatial_object@assays$SCT@scale.data=ScaleData(spatial_object,features = rownames(spatial_object))
  cat("saving RData file ... \n")
  save.image(paste0(filename,".RData"))
  return(spatial_object)
}





