## ST data analysis ----
source("E:/important/newest scripts/functions/my_spatial_analysis.R")
source("E:/important/newest scripts/functions/source.R")
source("E:/important/newest scripts/functions/spatial_findmarkers.R")
source("E:/important/newest scripts/functions/myspatial_featureplot.R")

### import the data ----
## p0928s1 ----
data.path="./20210928WF1"
p0928s1=my_spatial_analysis(data.path = data.path,colors = colors_crc_spatial,
                            filename = "p0928s1")
## p0928s2 ----
data.path="./20210928WF2"
p0928s2=my_spatial_analysis(data.path = data.path,colors = colors_crc_spatial,
                            filename = "p0928s2",resolution = 0.3,pt.size.factor = 0.8)
## p0928s3 ----
data.path="./20210928WF3-2"
p0928s3=my_spatial_analysis(data.path = data.path,colors = colors_crc_spatial,
                            filename = "p0928s3",resolution = 0.3)
## p1125s3b2 ----
data.path="./20211125WF3"
p1125s3b2=my_spatial_analysis(data.path = data.path,colors = colors_crc_spatial,
                              filename = "p1125s3b2",resolution = 0.2)
## p1125s4b2 ----
data.path="./20211125WF4_2"
p1125s4b2=my_spatial_analysis(data.path = data.path,colors = colors_crc_spatial,
                              filename = "p1125s4b2",resolution = 0.3)
## p1129s3 ----
data.path="./20211129WF3_3"
p1129s3=my_spatial_analysis(data.path = data.path,colors = colors_crc_spatial,
                            filename = "p1129s3",resolution = 0.3)


#### downstream analysis ----
## analysis for crc spatial ----
colors_sp=c("cyan","magenta","chartreuse","yellow4","gold","deepskyblue","purple","darkorange","deeppink")
## 1.p0928s1 for example ----
HE.p0928s1.path="./20210928WF1/spatial/tissue_lowres_image.png"
spotlight.sp.p0928s1=p0928s1
spotlight_decon_p0928s1 <- spotlight_deconvolution(se_sc = spotlight.object,
                                                   counts_spatial = spotlight.sp.p0928s1@assays$Spatial@counts,
                                                   clust_vr = "subclass",
                                                   cluster_markers = cluster_markers_all,
                                                   cl_n = 100,
                                                   hvg = 3000,
                                                   ntop = NULL,
                                                   transf = "uv",
                                                   method = "nsNMF",
                                                   min_cont = 0.05)
decon_mtrx_p0928s1 <- spotlight_decon_p0928s1[[2]]
spotlight.sp.p0928s1@meta.data <- cbind(spotlight.sp.p0928s1@meta.data, decon_mtrx_p0928s1)
## percent of clusters in regions ----
percentplot_spatial_decon(spotlight.sp.p0928s1,decon_matrix = decon_mtrx_p0928s1,
                          filename = "p0928s1_clusters",colors = colors_sp,height = 3)
## deconvolution pie plot ----
pdf("spotlight_p0928s1.pie.pdf",width = 12,height = 10)
SPOTlight::spatial_scatterpie(se_obj = spotlight.sp.p0928s1,
                              cell_types_all = cell_types_all,
                              img_path = HE.p0928s1.path,
                              pie_scale = 0.4)+scale_fill_manual(values = colors_sp)
dev.off()

## dominant cluster spatialplot ----
decon_mtrx_p0928s1=as.data.frame(decon_mtrx_p0928s1)
decon_mtrx_p0928s1$spotlight.cluster=colnames(decon_mtrx_p0928s1)[apply(decon_mtrx_p0928s1[,-ncol(decon_mtrx_p0928s1)],
                                                                        1,which.max)]
spotlight.sp.p0928s1$cluster=Idents(spotlight.sp.p0928s1)
spotlight.sp.p0928s1$spotlight.cluster=decon_mtrx_p0928s1$spotlight.cluster
spotlight.sp.p0928s1$spotlight.cluster=factor(spotlight.sp.p0928s1$spotlight.cluster,levels = cell_types_all)
Idents(spotlight.sp.p0928s1)="spotlight.cluster"
#Idents(spotlight.sp.p0928s1)="cluster"
pdf("p0928s1_spotlight_cluster.pdf",width = 6.5,height = 5)
SpatialDimPlot(spotlight.sp.p0928s1,pt.size.factor = 1.4)+scale_fill_manual(values = colors_sp)
dev.off()

p0928s1.label=spotlight.sp.p0928s1@meta.data[,c(16:23,26)]
write.csv(p0928s1.label,file = "p0928s1_label_transfer_bc.csv")
Idents(spotlight.sp.p0928s1)="cluster"

## spatial dimplot and umap plot ----
p1 <- DimPlot(spotlight.sp.p0928s1, reduction = "umap", label = F,cols = colors_crc)
p2 <- SpatialDimPlot(spotlight.sp.p0928s1,pt.size.factor = 1.4,
                     label = T,label.size = 5)+scale_fill_manual(values = colors_crc)
pdf(paste0("p0928s1","_spatialdimplot.pdf"),width = 6.5*1.8,height = 5)
p3=p1 + p2
print(p3)
dev.off()

## featureplot ----
myspatial_featureplot(spotlight.sp.p0928s1,filename = "p0928s1",genename = c("EPCAM","ALB"),
                      pt.size.factor = 1.4)

#############################################################
#########****** cpdb *****###################################
#############################################################
## p0928s1 for example ----
p0928s1_counts=as.data.frame(spotlight.sp.p0928s1@assays$SCT@data)
p0928s1_counts=data.frame(Gene=rownames(p0928s1_counts), p0928s1_counts)
colnames(p0928s1_counts)=str_replace_all(colnames(p0928s1_counts),pattern = "\\.","-")
p0928s1_meta=spotlight.sp.p0928s1@meta.data
p0928s1_meta=data.frame(Cell=rownames(p0928s1_meta),
                        cell_type=spotlight.sp.p0928s1$spotlight.cluster)
setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/p0928s1/cpdb")
write.table(p0928s1_counts, "p0928s1_counts.txt", row.names=F, sep='\t')
write.table(p0928s1_meta, "p0928s1_meta.txt", row.names=F, sep='\t')
## run cpdb, python code ----
#cellphonedb method statistical_analysis p0928s1_meta.txt p0928s1_counts.txt --counts-data=gene_name --iterations=100 --threads=64
setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/p0928s1/cpdb/out")
mypvals <- read.delim("pvalues.txt", check.names = FALSE)
mymeans <- read.delim("means.txt", check.names = FALSE)
sig.means=read.delim("significant_means.txt",check.names = F)
sig.means[is.na(sig.means)]=0

source("./cell_cluster_interaction.R")
mypvals=order_LR_pairs(mypvals)
mymeans=order_LR_pairs(mymeans)
sig.means=order_LR_pairs(sig.means,sig.matrix = T)

mypvals=mypvals[-which(duplicated(mypvals$interacting_pair)),]
mymeans=mymeans[-which(duplicated(mymeans$interacting_pair)),]
sig.means=sig.means[-which(duplicated(sig.means$interacting_pair)),]
rownames(sig.means)=sig.means$interacting_pair
sig.means=sig.means[mymeans$interacting_pair,]

write.table(mypvals,"pvalues_reorder.txt",sep = "\t",col.names = T)
write.table(mymeans,"means_reorder.txt",sep = "\t",col.names = T)
write.table(sig.means,"significant_means_reorder.txt",sep = "\t",col.names = T)

interaction_matrix=interaction_calculation(sig.matrix=sig.means)

bubble_matrix=interaction_bubble(pvalue.path = "pvalues_reorder.txt",means.path = "means_reorder.txt",filename = "p0928s1_LR_cpdb",
                                 pattern = "\\|",width = 20,
                                 height = 50,legend.midpoint = 1.5,means.cutoff = 0.5,pvalue.cutoff = Inf,
                                 vjust = 0.5,
                                 hjust = 0)
pattern="Tumor|Fibroblast"
choose_bubble_matrix=interaction_bubble_choose(pvalue.path = "pvalues_reorder.txt",means.path = "means_reorder.txt",
                                               filename = "p0928s1_choose_LR",
                                               pattern = pattern,width = 10,height = 10,
                                               choose.interaction = choose_interaction,
                                               legend.midpoint = 1.5,pvalue.cutoff = Inf,
                                               means.cutoff = 1.5,scale.size = 4,vjust = 0.5,
                                               hjust = 0)


######################################################
#######******* featureplot *******####################
######################################################
setwd("./featureplot")
png(file = paste0("p0928s1_featureplot.png"),width = 3000,height = 500)
SpatialFeaturePlot(spotlight.sp.p0928s1,
                   features = c("EPCAM","ALB","COL1A1","CD79A","CD3D","CD14"),ncol = 6)
dev.off()

png(file = paste0("p0928s2_featureplot.png"),width = 3000,height = 500)
SpatialFeaturePlot(spotlight.sp.p0928s2,
                   features = c("EPCAM","ALB","COL1A1","CD79A","CD3D","CD14"),ncol = 6)
dev.off()

png(file = paste0("p0928s3_featureplot.png"),width = 3000,height = 500)
SpatialFeaturePlot(spotlight.sp.p0928s3,
                   features = c("EPCAM","ALB","COL1A1","CD79A","CD3D","CD14"),ncol = 6)
dev.off()

png(file = paste0("p1125s3b2_featureplot.png"),width = 3000,height = 500)
SpatialFeaturePlot(spotlight.sp.p1125s3b2,
                   features = c("EPCAM","ALB","COL1A1","CD79A","CD3D","CD14"),ncol = 6)
dev.off()

png(file = paste0("p1125s4b2_featureplot.png"),width = 3000,height = 500)
SpatialFeaturePlot(spotlight.sp.p1125s4b2,
                   features = c("EPCAM","ALB","COL1A1","CD79A","CD3D","CD14"),ncol = 6)
dev.off()

png(file = paste0("p1129s3_featureplot.png"),width = 3000,height = 500)
SpatialFeaturePlot(spotlight.sp.p1129s3,
                   features = c("EPCAM","ALB","COL1A1","CD79A","CD3D","CD14"),ncol = 6)
dev.off()

######################################################
#######******* tumor region *******###################
######################################################
setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/tumor_region")
spotlight.sp.p0928s1$tumor_region=ifelse(spotlight.sp.p0928s1$cluster%in%c("0","7"),"PT","T")
Idents(spotlight.sp.p0928s1)="tumor_region"
png(file = paste0("p0928s1_tumor_region.png"),width = 1600,height = 1200,res=600)
SpatialDimPlot(spotlight.sp.p0928s1,image.alpha = 0,stroke = 0,pt.size.factor = 2)+scale_fill_manual(values = colors2)
dev.off()

spotlight.sp.p0928s2$tumor_region=ifelse(spotlight.sp.p0928s2$cluster%in%c("0","3","5","6"),"T","PT")
Idents(spotlight.sp.p0928s2)="tumor_region"
png(file = paste0("p0928s2_tumor_region.png"),width = 1600,height = 1200,res=600)
SpatialDimPlot(spotlight.sp.p0928s2,image.alpha = 0,stroke = 0,pt.size.factor = 1.7)+scale_fill_manual(values = colors2[c(2,1)])
dev.off()

spotlight.sp.p0928s3$tumor_region=ifelse(spotlight.sp.p0928s3$cluster%in%c("3","4"),"PT","T")
Idents(spotlight.sp.p0928s3)="tumor_region"
png(file = paste0("p0928s3_tumor_region.png"),width = 1600,height = 1200,res=600)
SpatialDimPlot(spotlight.sp.p0928s3,image.alpha = 0,stroke = 0,pt.size.factor = 1.8)+scale_fill_manual(values = colors2[c(2,1)])
dev.off()

spotlight.sp.p1125s3b2$tumor_region=ifelse(spotlight.sp.p1125s3b2$cluster%in%c("0","2","5"),"T","PT")
Idents(spotlight.sp.p1125s3b2)="tumor_region"
png(file = paste0("p1125s3b2_tumor_region.png"),width = 1600,height = 1200,res=600)
SpatialDimPlot(spotlight.sp.p1125s3b2,image.alpha = 0,stroke = 0,pt.size.factor = 1.6)+scale_fill_manual(values = colors2[c(2,1)])
dev.off()

spotlight.sp.p1125s4b2$tumor_region=ifelse(spotlight.sp.p1125s4b2$cluster%in%c("1","4"),"PT","T")
Idents(spotlight.sp.p1125s4b2)="tumor_region"
png(file = paste0("p1125s4b2_tumor_region.png"),width = 1600,height = 1200,res=600)
SpatialDimPlot(spotlight.sp.p1125s4b2,image.alpha = 0,stroke = 0,pt.size.factor = 1.7)+scale_fill_manual(values = colors2)
dev.off()

spotlight.sp.p1129s3$tumor_region=ifelse(spotlight.sp.p1129s3$cluster%in%c("4"),"T","PT")
Idents(spotlight.sp.p1129s3)="tumor_region"
png(file = paste0("p1129s3_tumor_region.png"),width = 1600,height = 1200,res=600)
SpatialDimPlot(spotlight.sp.p1129s3,image.alpha = 0,stroke = 0,pt.size.factor = 1.7)+scale_fill_manual(values = colors2[c(2,1)])
dev.off()

### GSVA of hallmark genesets ----
gspath="/mnt/ssd1/longjie/newest_scripts/database/gsea/h.all.v7.1.symbols.gmt"
geneset=getGmt(gspath)

em.p0928s1=AverageExpression(spotlight.sp.p0928s1)$SCT
em.p0928s2=AverageExpression(spotlight.sp.p0928s2)$SCT
em.p0928s3=AverageExpression(spotlight.sp.p0928s3)$SCT
em.p1125s3b2=AverageExpression(spotlight.sp.p1125s3b2)$SCT
em.p1125s4b2=AverageExpression(spotlight.sp.p1125s4b2)$SCT
em.p1129s3=AverageExpression(spotlight.sp.p1129s3)$SCT

intergenes=intersect(rownames(em.p0928s1),intersect(rownames(em.p0928s2),intersect(rownames(em.p0928s3),
                                                                                   intersect(rownames(em.p1125s3b2),intersect(rownames(em.p1125s4b2),
                                                                                                                              rownames(em.p1129s3))))))
em.merge=cbind(em.p0928s1[intergenes,])
setdiff(rownames(em.p0928s1),intergenes)
escore.p0928s1 <- gsva(as.matrix(em.p0928s1), min.sz=2, max.sz=2000,
                       geneset,verbose=FALSE, parallel.sz=8,method="ssgsea")
escore.p0928s1=as.data.frame(escore.p0928s1)
colnames(escore.p0928s1)=paste0("C1_",colnames(escore.p0928s1))

escore.p0928s2 <- gsva(as.matrix(em.p0928s2), min.sz=2, max.sz=2000,
                       geneset,verbose=FALSE, parallel.sz=8,method="ssgsea")
escore.p0928s2=as.data.frame(escore.p0928s2)
colnames(escore.p0928s2)=paste0("C2_",colnames(escore.p0928s2))

escore.p0928s3 <- gsva(as.matrix(em.p0928s3), min.sz=2, max.sz=2000,
                       geneset,verbose=FALSE, parallel.sz=8,method="ssgsea")
escore.p0928s3=as.data.frame(escore.p0928s3)
colnames(escore.p0928s3)=paste0("C3_",colnames(escore.p0928s3))

escore.p1125s3b2 <- gsva(as.matrix(em.p1125s3b2), min.sz=2, max.sz=2000,
                         geneset,verbose=FALSE, parallel.sz=8,method="ssgsea")
escore.p1125s3b2=as.data.frame(escore.p1125s3b2)
colnames(escore.p1125s3b2)=paste0("L1_",colnames(escore.p1125s3b2))

escore.p1125s4b2 <- gsva(as.matrix(em.p1125s4b2), min.sz=2, max.sz=2000,
                         geneset,verbose=FALSE, parallel.sz=8,method="ssgsea")
escore.p1125s4b2=as.data.frame(escore.p1125s4b2)
colnames(escore.p1125s4b2)=paste0("L2_",colnames(escore.p1125s4b2))

escore.p1129s3 <- gsva(as.matrix(em.p1129s3), min.sz=2, max.sz=2000,
                       geneset,verbose=FALSE, parallel.sz=8,method="ssgsea")
escore.p1129s3=as.data.frame(escore.p1129s3)
colnames(escore.p1129s3)=paste0("C4_",colnames(escore.p1129s3))

escore.merge=cbind(escore.p0928s1,cbind(escore.p0928s2,cbind(escore.p0928s3,cbind(escore.p1129s3,
                                                                                  cbind(escore.p1125s3b2,escore.p1125s4b2)))))
escore_scale=as.data.frame(t(scale(as.data.frame(t(escore.merge)))))

rownames(escore_scale)=str_replace_all(rownames(escore_scale),"HALLMARK_","")
rownames(escore_scale)=str_to_title(rownames(escore_scale))

annotation_col=data.frame(Patients=str_split_fixed(colnames(escore_scale),"_",2)[,1],
                          Module=ifelse(colnames(escore_scale)%in%c("L2_T","L2_PT"),
                                        "Module1",
                                        ifelse(colnames(escore_scale)%in%c("L1_PT",
                                                                           "C2_PT","C4_PT"),
                                               "Module3","Module2")))
rownames(annotation_col)=colnames(escore_scale)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
png(paste0("GSVA_heatmap_regeion_hallmarker.png"),
    width = 3600,height = 4000,res=600)
p=pheatmap(escore_scale,cluster_rows = T,cluster_cols = T,
           treeheight_row = 0,fontsize_row = 8,fontsize_col = 8,
           color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           border_color = "white",breaks = bk,annotation_col = annotation_col)
print(p)
dev.off()

###################################################
###******       featureplot figures      *****#####
###################################################
p1=SpatialFeaturePlot(spotlight.sp.p0928s1,features = "F3")
p2=SpatialFeaturePlot(spotlight.sp.p0928s2,features = "F3")
p3=SpatialFeaturePlot(spotlight.sp.p0928s3,features = "F3")
p4=SpatialFeaturePlot(spotlight.sp.p1125s3b2,features = "F3")
p5=SpatialFeaturePlot(spotlight.sp.p1125s4b2,features = "F3")
p6=SpatialFeaturePlot(spotlight.sp.p1129s3,features = "F3")

fix.sc <- scale_fill_gradientn(colors = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick")
                               ,limits = c(0.0,3.5),breaks = c(0.0,1.0,2.0,3.0))
p <- lapply(list(p1,p2,p3,p6,p4,p5), function (x) x + fix.sc)
setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/Featureplot_F3")
png(paste0("Featureplot_F3.png"),
    width = 3600,height = 2500,res=600)
CombinePlots(p)
dev.off()

p1=SpatialFeaturePlot(spotlight.sp.p0928s1,features = "NRG1")
p2=SpatialFeaturePlot(spotlight.sp.p0928s2,features = "NRG1")
p3=SpatialFeaturePlot(spotlight.sp.p0928s3,features = "NRG1")
p4=SpatialFeaturePlot(spotlight.sp.p1125s3b2,features = "NRG1")
p5=SpatialFeaturePlot(spotlight.sp.p1125s4b2,features = "NRG1")
p6=SpatialFeaturePlot(spotlight.sp.p1129s3,features = "NRG1")

fix.sc <- scale_fill_gradientn(colors = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick")
                               ,limits = c(0.0,3.5),breaks = c(0.0,1.0,2.0,3.0))
p <- lapply(list(p1,p2,p3,p6,p4,p5), function (x) x + fix.sc)
setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/Featureplot_F3")
png(paste0("Featureplot_NRG1.png"),
    width = 3600,height = 2500,res=600)
CombinePlots(p)
dev.off()

p1=SpatialFeaturePlot(spotlight.sp.p0928s1,features = "ERBB3")
p2=SpatialFeaturePlot(spotlight.sp.p0928s2,features = "ERBB3")
p3=SpatialFeaturePlot(spotlight.sp.p0928s3,features = "ERBB3")
p4=SpatialFeaturePlot(spotlight.sp.p1125s3b2,features = "ERBB3")
p5=SpatialFeaturePlot(spotlight.sp.p1125s4b2,features = "ERBB3")
p6=SpatialFeaturePlot(spotlight.sp.p1129s3,features = "ERBB3")
CombinePlots(list(p1,p2,p3,p4,p5,p6))

fix.sc <- scale_fill_gradientn(colors = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick")
                               ,limits = c(0.0,3.5),breaks = c(0.0,1.0,2.0,3.0))
p <- lapply(list(p1,p2,p3,p6,p4,p5), function (x) x + fix.sc)
setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/Featureplot_F3")
png(paste0("Featureplot_ERBB3.png"),
    width = 3600,height = 2500,res=600)
CombinePlots(p)
dev.off()

## correlation of genes-----
p0928s1.matrix=as.data.frame(t(spotlight.sp.p0928s1@assays$Spatial@data))
p1=ggplot(p0928s1.matrix,aes(F3,NRG1))+geom_point(color="slateblue",size=2)+
  theme(panel.border = element_rect(fill = "NA",size = 1,color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank())+
  geom_smooth(method= "lm",color="firebrick2")+
  stat_cor(method = "pearson")+labs(title = "C1_correlation")+
  theme(plot.title = element_text(hjust = 0.5))

png(paste0("Correlation_F3_NRG1.png"),
    width = 7800,height = 5000,res=600)
print(p1)
dev.off()


## correlation of scores --
F02_markers=FindMarkers(tumor.combined,ident.1 = "F02_fibrblast_MCAM",logfc.threshold = 2,
                        min.pct = 0.5,only.pos = T)
F02_markers=F02_markers[order(F02_markers$avg_log2FC,decreasing = T),]
F02_markers_genename=list(rownames(F02_markers))

CD8_markers=FindMarkers(immune.combined,ident.1 = "T03_CD8_CXCL13",logfc.threshold = 1,
                        min.pct = 0.25,only.pos = T)
CD8_markers=CD8_markers[order(CD8_markers$avg_log2FC,decreasing = T),]
T03_markers_genename=list(rownames(CD8_markers))

spotlight.sp.p1125s3b2=AddModuleScore(spotlight.sp.p1125s3b2,features = F02_markers_genename,name = "F02_MCAM")
spotlight.sp.p1125s4b2=AddModuleScore(spotlight.sp.p1125s4b2,features = F02_markers_genename,name = "F02_MCAM")

spotlight.sp.p1125s3b2=AddModuleScore(spotlight.sp.p1125s3b2,features = T03_markers_genename,name = "CD8_CXCL13")
spotlight.sp.p1125s4b2=AddModuleScore(spotlight.sp.p1125s4b2,features = T03_markers_genename,name = "CD8_CXCL13")

p1=SpatialFeaturePlot(spotlight.sp.p1125s3b2,features = "F02_MCAM1",min.cutoff = "q5",pt.size.factor = 1.4)
p2=SpatialFeaturePlot(spotlight.sp.p1125s3b2,features = "CD8_CXCL131",min.cutoff = "q5",pt.size.factor = 1.4)

p3=SpatialFeaturePlot(spotlight.sp.p1125s4b2,features = "F02_MCAM1",min.cutoff = "q5",pt.size.factor = 1.4)
p4=SpatialFeaturePlot(spotlight.sp.p1125s4b2,features = "CD8_CXCL131",min.cutoff = "q5",pt.size.factor = 1.4)
setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/Featureplot_F3")
png(paste0("Featureplot_F02_MCAM_CD8_CXCL13.png"),
    width = 4000,height = 3200,res=600)
CombinePlots(list(p1,p2,p3,p4))
dev.off()

spotlight.sp.p1125s3b2.matrix=rbind(F02_MCAM=spotlight.sp.p1125s3b2$F02_MCAM1,
                                    CD8_CXCL13=spotlight.sp.p1125s3b2$CD8_CXCL131)
spotlight.sp.p1125s3b2.matrix=as.data.frame(t(spotlight.sp.p1125s3b2.matrix))

spotlight.sp.p1125s4b2.matrix=rbind(F02_MCAM=spotlight.sp.p1125s4b2$F02_MCAM1,
                                    CD8_CXCL13=spotlight.sp.p1125s4b2$CD8_CXCL131)
spotlight.sp.p1125s4b2.matrix=as.data.frame(t(spotlight.sp.p1125s4b2.matrix))


p1=ggplot(spotlight.sp.p1125s3b2.matrix,aes(F02_MCAM,CD8_CXCL13))+geom_point(color="slateblue",size=2)+
  theme(panel.border = element_rect(fill = "NA",size = 1,color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank())+
  geom_smooth(method= "lm",color="firebrick2")+
  #stat_cor(method = "pearson")+labs(title = "L1_correlation")+
  theme(plot.title = element_text(hjust = 0.5))

p2=ggplot(spotlight.sp.p1125s4b2.matrix,aes(F02_MCAM,CD8_CXCL13))+geom_point(color="slateblue",size=2)+
  theme(panel.border = element_rect(fill = "NA",size = 1,color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank())+
  geom_smooth(method= "lm",color="firebrick2")+
  #stat_cor(method = "pearson")+labs(title = "L2_correlation")+
  theme(plot.title = element_text(hjust = 0.5))

setwd("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/Featureplot_F3")
png(paste0("Correlation_F02_MCAM_CD8_CXCL13.png"),
    width = 6500,height = 3200,res=600)
CombinePlots(list(p1,p2))
dev.off()

png(paste0("Correlation_F02_MCAM_CD8_CXCL13_nopvalue.png"),
    width = 6500,height = 3200,res=600)
CombinePlots(list(p1,p2))
dev.off()

### GSVA pheatmap for tumor and normal ----
gspath="/mnt/ssd1/longjie/newest_scripts/database/gsea/h.all.v7.1.symbols.gmt"
geneset=getGmt(gspath)

em.p0928s1=AverageExpression(spotlight.sp.p0928s1)$SCT
em.p0928s2=AverageExpression(spotlight.sp.p0928s2)$SCT
em.p0928s3=AverageExpression(spotlight.sp.p0928s3)$SCT
em.p1125s3b2=AverageExpression(spotlight.sp.p1125s3b2)$SCT
em.p1125s4b2=AverageExpression(spotlight.sp.p1125s4b2)$SCT
em.p1129s3=AverageExpression(spotlight.sp.p1129s3)$SCT

intergenes=intersect(rownames(em.p0928s1),intersect(rownames(em.p0928s2),intersect(rownames(em.p0928s3),
                                                                                   intersect(rownames(em.p1125s3b2),intersect(rownames(em.p1125s4b2),
                                                                                                                              rownames(em.p1129s3))))))
colnames(em.p0928s1)=paste0("C1_",colnames(em.p0928s1))
colnames(em.p0928s2)=paste0("C2_",colnames(em.p0928s2))
colnames(em.p0928s3)=paste0("C3_",colnames(em.p0928s3))
colnames(em.p1125s3b2)=paste0("L1_",colnames(em.p1125s3b2))
colnames(em.p1125s4b2)=paste0("L2_",colnames(em.p1125s4b2))
colnames(em.p1129s3)=paste0("C4_",colnames(em.p1129s3))

em.merge=cbind(em.p0928s1[intergenes,],cbind(em.p0928s2[intergenes,],cbind(em.p0928s3[intergenes,],
                                                                           cbind(em.p1125s3b2[intergenes,],
                                                                                 cbind(em.p1125s4b2[intergenes,],
                                                                                       em.p1129s3[intergenes,])))))

escore.merge <- gsva(as.matrix(em.merge), min.sz=2, max.sz=2000,
                     geneset,verbose=FALSE, parallel.sz=8,method="ssgsea")

escore_scale=as.data.frame(t(scale(as.data.frame(t(escore.merge)))))

rownames(escore_scale)=str_replace_all(rownames(escore_scale),"HALLMARK_","")
rownames(escore_scale)=str_to_title(rownames(escore_scale))
pheatmap(escore_scale)

annotation_col=data.frame(Patients=str_split_fixed(colnames(escore_scale),"_",2)[,1],
                          Module=ifelse(colnames(escore_scale)%in%c("L2_T","L2_PT"),
                                        "Module1",
                                        ifelse(colnames(escore_scale)%in%c("L1_PT",
                                                                           "C2_PT","C4_PT"),
                                               "Module3","Module2")))
rownames(annotation_col)=colnames(escore_scale)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
png(paste0("GSVA_heatmap_regeion_hallmarker.png"),
    width = 3600,height = 4000,res=600)
p=pheatmap(escore_scale,cluster_rows = T,cluster_cols = T,
           treeheight_row = 0,fontsize_row = 8,fontsize_col = 8,
           color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           border_color = "white",breaks = bk,annotation_col = annotation_col)
print(p)
dev.off()
getwd()


save.image("/mnt/ssd1/longjie/help_analysis/wangfei_SA/20230205_spatial_crc/20230211_CRC_spatial.RData")











