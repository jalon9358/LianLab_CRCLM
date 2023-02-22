## CRC figures ----
## Figure 1B -----
pdf("Figure 1B.allnonimmune_dimplot.pdf",width = width.ppi*1.5,height = height.ppi)
DimPlot(object = tumor.combined,reduction = 'umap',label = F,cols = colors_tumor)
dev.off()

## Figure 1C ----
pdf("Figure 1C.small_tumor_cluster_diff.pdf",width = width.ppi*1.5,height = height.ppi*1.5)
print(volcano_tumor_anno)
dev.off()

pdf("Figure 1C.small_tumor_cluster_diff_nolabel.pdf",width = width.ppi*1.2,height = height.ppi*1.2)
print(volcano_tumor)
dev.off()

## Figure S1A ----
DefaultAssay(tumor.combined)="RNA"
featureplot_rna=c("EPCAM","COL1A1","PECAM1","SOX9","COL1A2","CD34")
myfeatureplot(tumor.combined,filename = "Figure S1A.tumor_RNA",genename = featureplot_rna,reduction = "umap",
              cols = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick"))

## Figure S1C ----
filename="alltumor_small"
pdf(paste0("Figure S1C.","heatmap_",filename,"_width_top5.pdf"),width =13, height = 10)
f1=DoHeatmap(tumor_ave_cluster, features = heatgene,
             label = TRUE,size=3,assay = "RNA",lines.width = 1,
             group.colors = colors_tumor,draw.lines = F)+
  scale_fill_gradientn(colors =  rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
print(f1)
dev.off()

## Figure S1E ----
percentplot_choose(seurat_object = tumor.combined,filename = "Figure S1E.small_tumor_organ_patients",
                   choose = "patients_organ",width = 10,height = 4,legend.key.size = 4,
                   angle = 90,colors = colors_tumor)

## Figure S2A ----
pdf("Figure S2A.tumorcell_dimplot.pdf",width = width.ppi*1.1,height = height.ppi)
DimPlot(object = tumorcell.combined, reduction = 'umap',label = F,cols = colors_tumorcell)
dev.off()

pdf("Figure S2A.tumorcell_organs_dimplot.pdf",width = width.ppi*1.8,height = height.ppi)
DimPlot(object = tumorcell.combined, reduction = 'umap',label = F,cols = colors_tumorcell,raster = F,
        split.by = "organs",ncol = 2)
dev.off()
FeaturePlot(immune.combined,features = "MS4A2")

## Figure S2C ----
levels(fibroblast.combined)=sort(levels(fibroblast.combined))
setwd("/GPUFS/scut_zxlian_4/longjie/CRC_20220523_fibroblast")
pdf("Figure S2D.fibroblast_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = fibroblast.combined, reduction = 'umap',label = F,cols = colors_fibroblast)
dev.off()

pdf("Figure S2D.fibroblast_dimplot_organs.pdf",width = width.ppi*1.8,height = height.ppi*1)
DimPlot(object = fibroblast.combined, reduction = 'umap',label = F,cols=colors_fibroblast,
        split.by = "organs",
        ncol = 2)
dev.off()
FeaturePlot(immune.combined,features = "CXCR1")
FeaturePlot(fibroblast.combined,features = "MCAM",split.by = "organs")
FeaturePlot(tumorcell.combined,features = "F3",split.by = "organs")

dev.off()
## Figure S2D ----
pdf("Figure S2D.fibroblast_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = fibroblast.combined, reduction = 'umap',label = F,cols = colors_fibroblast)
dev.off()

## Figure S2F ----
pdf("Figure S2F.fibroblast_dimplot_patients.pdf",width = width.ppi*3.8,height = height.ppi*2)
DimPlot(object = fibroblast.combined, reduction = 'umap',label = F,
        cols = colors_fibroblast,pt.size = 1.5)+
  facet_wrap(fibroblast.combined$patients~fibroblast.combined$organs,ncol = 5)
dev.off()

## Figure 2A ----
percentplot_choose(seurat_object = tumorcell.combined,
                   filename = "Figure 2A.tumorcell_organ_patients",choose = "patients_organ",
                   width = 11.2,height = 4,legend.key.size = 4,angle = 90,colors = colors_tumorcell)
percentplot_choose(seurat_object = tumorcell.combined,
                   filename = "Figure 2A.tumor_organs",choose = "organs",
                   width = 4,height = 4,legend.key.size = 4,angle = 90,colors = colors_tumorcell)

## Figure 2B ----
pdf(paste0("Figure 2B.GSVA_cluster_diff_","crc_tumorcell",".pdf"),width = 1*6.5,height = 6)
p=pheatmap(tumorcell_gsva,fontsize_row = 6,border_color = "black",
           cluster_rows = F,cluster_cols = F,treeheight_row = 0,
           color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
print(p)
dev.off()

pdf("Figure 2B.tumorcell_transcription_net_scenic.pdf",width = width.ppi*1.1,height = height.ppi*2)
p=pheatmap(regulonActivity_byCellType_Scaled,border_color = "black",
           cluster_cols = F,treeheight_row = 0,
           color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
print(p)
dev.off()

## Figure 2C ----
stacked_violin_plot(gene = c("LGR5"),col = colors_tumorcell,
                    filename = "Figure 2C.LGR5",seurat_object = tumorcell.combined,
                    Mean = T,flip = F,width = 10,height = 4,legend.position = "right",
                    limits.max = 1.5)
stacked_violin_plot(gene = c("EPCAM"),col = colors_tumorcell,
                    filename = "Figure 2C.EPCAM",seurat_object = tumorcell.combined,
                    Mean = T,flip = F,width = 10,height = 4,legend.position = "right",
                    limits.max = 7)
stacked_violin_plot(gene = c("CDH1"),col = colors_tumorcell,
                    filename = "Figure 2C.CDH1",seurat_object = tumorcell.combined,
                    Mean = T,flip = F,width = 10,height = 4,legend.position = "right",
                    limits.max = 3)
## Figure 2D ----
percentplot_choose(seurat_object = fibroblast.combined,
                   filename = "Figure 2D.fibroblast_organ_patients",choose = "patients_organ",
                   width = 11.2,height = 4,legend.key.size = 4,angle = 90,colors = colors_fibroblast)
percentplot_choose(seurat_object = fibroblast.combined,
                   filename = "Figure 2D.fibroblast_organs",choose = "organs",
                   width = 4,height = 4,legend.key.size = 4,angle = 90,colors = colors_fibroblast)
FeaturePlot(tumor.combined,features = "MMP2")
## Figure 3G ----
myfeatureplot(myeloid.combined,genename =c("LAMP3","CCR7","CD40","CD80","CD86","HLA.DRA"),
              filename = "Figure 3G.myeloid",reduction = "umap",
              cols = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick"))
myfeatureplot(tumor.combined,genename =c("MMP2"),
              filename = "MMP2",reduction = "umap",
              cols = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick"))
## Figure 3H ----
stacked_violin_plot(gene = c("CD40","CD80","CD86"),col = colors_myeloid,
                    filename = "Figure 3H.CD40",seurat_object = myeloid.combined,
                    Mean = T,flip = F,width = 10,height = 4,legend.position = "right",
                    limits.max = 3)

## Figure S3B ----
myeloid_markers=c("CD14","CXCL9","CD1C","ZNF331","SPP1","HSPA6",
                  "SERPINB2","MKI67","CPNE3","FCGR3A","LAMP3")
myfeatureplot(myeloid.combined,genename =myeloid_markers,
              filename = "Figure S3B.myeloid",reduction = "umap",
              cols = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick"))

## Figure 4C ----
pdf("Figure 4C.T_cluster_diff_nolabel.pdf",width = 6.5,height = 5)
print(fig_diff_t_nolabel)
dev.off()

pdf("Figure 4C.T_cluster_diff.pdf",width = 6.5,height = 5)
print(fig_diff_t_label)
dev.off()

## Figure 4F ----
pdf(paste0("Figure 4F.GSVA_cluster_diff_","crc_CD8",".pdf"),width = 1*5.5,height = 6)
p=pheatmap(cd8_gsva,fontsize_row = 6,cluster_rows = F,cluster_cols = F,
           treeheight_row = 0,border_color = "black",
           color = colorRampPalette(colors = c("blue4","white","red"))(100))
print(p)
dev.off()

## Figure 4G ----
pdf(paste0("Figure 4G.GSVA_cluster_diff_","crc_CD4",".pdf"),width = 1*6,height = 6)
p=pheatmap(cd4_gsva,fontsize_row = 6,cluster_rows = F,cluster_cols = F,
           treeheight_row = 0,border_color = "black",
           color = colorRampPalette(colors = c("blue4","white","red"))(100))
print(p)
dev.off()

## Figure S4B ----
pdf("Figure S4B.CD4_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = cd4.combined, reduction = 'umap',label = F,cols = colors_cd4)
dev.off()

pdf("Figure S4B.CD4_dimplot_samples.pdf",width = width.ppi*1.1,height = height.ppi*1.5)
DimPlot(object = cd4.combined, reduction = 'umap',label = F,cols = colors_cd4,split.by = "organs",
        ncol = 2)
dev.off()

## Figure S4C ----
filename="allcd8_small"
pdf(paste0("Figure S4C.heatmap_",filename,"_top5.pdf"),width =9, height = 8)
f1=DoHeatmap(cd8_ave_cluster, features = heatgene,
             label = TRUE,size=3,assay = "RNA",lines.width = 1,
             group.colors = colors_cd8,draw.lines = F)+
  scale_fill_gradientn(colors =  rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
print(f1)
dev.off()

## Figure S4D ----
filename="allcd4_small"
pdf(paste0("Figure S4D.heatmap_",filename,"_top5.pdf"),width =9, height = 8)
f1=DoHeatmap(cd4_ave_cluster, features = heatgene_cd4,
             label = TRUE,size=3,assay = "RNA",lines.width = 1,
             group.colors = colors_cd4,draw.lines = F)+
  scale_fill_gradientn(colors =  rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
print(f1)
dev.off()

## Figure S4E ----
DefaultAssay(cd8.combined)="antibody"
myfeatureplot(cd8.combined,genename ="CD69.pAbO",filename = "Figure S4E.CD8_CD69",
              reduction = "umap",cols = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick"),
              pt.size = 0.5)
DefaultAssay(cd8.combined)="RNA"

## Figure 5A ----
myvolcano(CD8_CXCL13vsOthers,gene.plot = 150,filename = "Figure 5A.CD8_diffgene_labeled",
          logfc.cutoff = 0.6,text.repel = T,width = 26,height = 15,text.size = 1)
myvolcano(CD8_CXCL13vsOthers,gene.plot = 150,filename = "Figure 5A.CD8_diffgene_nolabel",
          logfc.cutoff = 0.6,text.repel = F,pt.size = 2)
myvolcano(CD8_CXCL13vsOthers,gene.plot = 150,filename = "Figure 5A.CD8_diffgene_nolabel_nogrid",
          logfc.cutoff = 0.6,text.repel = F,pt.size = 2,grid = F)

## Figure 5B ----
DefaultAssay(cd8.combined)="antibody"
feature_gene=c("CD279.pAbO","TIM3.pAbO","LAG3.pAbO")
myfeatureplot(cd8.combined,genename =feature_gene,filename = "Figure 5B.CD8_function_antibody",
              reduction = "umap",cols = c("slateblue4","lightgreen","lightgoldenrod1",
                                          "darkorange","firebrick"))
DefaultAssay(cd8.combined)="RNA"

## Figure 5C ----
load("Figure 5C.CD8_choose_genecds.RData")
filename = "Figure 5C.CD8_choose_gene"
p1=plot_cell_trajectory(export$cds,color_by="celltype", cell_size=0.5,show_backbone=TRUE)+
  scale_color_manual(values = colors_cd8)
ggsave(paste0(filename,"monocle.celltype.pdf"),p1,width = width.ppi,height = height.ppi)

export$cds <- orderCells(export$cds, root_state = 9) 
p2=plot_cell_trajectory(export$cds,color_by="Pseudotime", cell_size=0.5,show_backbone=TRUE) 
ggsave(paste0(filename,"monocle.pseudotime.pdf"),p2,width = width.ppi,height = height.ppi)

## Figure 6C ----
pdf(paste0("Figure 6C.dotplot of NOTCH gene.pdf"),width = 1*6.5,height = 5)
DotPlot(t.combined,features = notch_gene,cols = c("orange3","darkviolet"))
dev.off()

## Figure 6D ----
pdf(paste0("Figure 6D.dotplot of NOTCH ligand gene.pdf"),width = 1*6.5,height = 5)
DotPlot(tumor.combined,features = notch_ligandgene,cols = c("orange3","darkviolet"))
dev.off()
levels(tumor.combined)=rev(levels(tumor.combined))

## Figure S5A ----
feature_gene=c("DLL1","DLL3","DLL4","JAG1","JAG2")
myfeatureplot(immune.combined,filename = "Figure S5A.immunecell_NOTCH_ligand",
              genename = feature_gene,reduction = "umap",
              cols = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick"))

myfeatureplot(tumor.combined,filename = "Figure S5A.tumor_NOTCH_ligand",
              genename = feature_gene,reduction = "umap",
              cols = c("slateblue4","lightgreen","lightgoldenrod1","darkorange","firebrick"))

## Figure 6E ----
source("/GPUFS/scut_zxlian_4/longjie/newest scripts/functions/cell_cluster_interaction.R")
notch_interaction=grep("DLL[1|3|4]_NOTCH|JAG[1|2]_NOTCH",mypvals$interacting_pair,value = T)
notch_interaction_DLL=grep("DLL[1|3|4]_NOTCH",mypvals$interacting_pair,value = T)
notch_bubble_matrix=interaction_bubble_choose(pvalue.path = "pvalues.txt",means.path = "means.txt",
                                              filename = "Figure 6E.NOTCH_cxcl13_cd8",
                                              pattern = "^Tu.*\\|T.*_CD8_CXCL13.*|^E.*\\|T.*_CD8_CXCL13.*|^F.*\\|T.*_CD8_CXCL13.*",
                                              choose.interaction = notch_interaction,
                                              legend.midpoint = 0.3,pvalue.cutoff = 10,
                                              means.cutoff = 0.2,scale.size = 10)
## Figure S5B ----
notch_cd4_bubble_matrix=interaction_bubble_choose(pvalue.path = "pvalues.txt",means.path = "means.txt",
                                                  filename = "Figure S5B.NOTCH_cxcl13_cd4",
                                                  pattern = "^Tu.*\\|T.*_CD4_CXCL13.*|^E.*\\|T.*_CD4_CXCL13.*|^F.*\\|T.*_CD4_CXCL13.*",
                                                  choose.interaction = notch_interaction,
                                                  legend.midpoint = 0.3,pvalue.cutoff = 10,
                                                  means.cutoff = 0.2,scale.size = 10)
immune.combined@assays$RNA@data[1:5,1:5]

## Figure S5G ----
source("/GPUFS/scut_zxlian_4/longjie/newest scripts/functions/mycellchat.R")
source("/GPUFS/scut_zxlian_4/longjie/newest scripts/functions/myrankNet.R")
#color.use=c("darkmagenta","brown1")
color.use=c("darkgoldenrod","darkmagenta")
gg1 <- myrankNet(merge.cellchat, mode = "comparison", stacked = T, do.stat = TRUE,
                 color.use = color.use,myp.cutoff = 1e-15)
gg2 <- myrankNet(merge.cellchat, mode = "comparison", stacked = F, do.stat = TRUE,
                 color.use = color.use,myp.cutoff = 1e-15)

pdf("Figure S5G.funtional_model_diff_p2.pdf",width = 8,height = 6)
gg1 + gg2
dev.off()

## Figure 6B ----
setwd("/GPUFS/scut_zxlian_4/longjie/CRC_20220613")
fibro_topgene=myfindmarkers(fibroblast.combined,filename = "Figure 6B.fibro_marker",colors = colors_fibroblast)
source("/GPUFS/scut_zxlian_4/longjie/newest scripts/functions/myviolin_sina.R")
myviolin_sina(fibroblast.combined,colors=colors_fibroblast,violin_gene = "NRG1")

## Figure 6. fibroblast pyscenic ----
set.seed(15454)
cellname=sample(colnames(fibroblast.combined),1500)
setwd("/GPUFS/scut_zxlian_4/longjie/CRC_20220613/pyscenic")
fib.scenic.choose=subset(fibroblast.combined,cells=cellname)
write.csv(t(as.matrix(fib.scenic.choose@assays$RNA@counts)), file = "fibroblast_choose.csv")

## Figure 7. volcano fibroblast ----
F_MCAMvsF_F3=FindMarkers(fibroblast.combined,ident.1 = "F02_fibrblast_MCAM",
                         ident.2 = "F04_fibroblast_F3")
F_MCAMvsF_F3$p_val_adj[which(F_MCAMvsF_F3$p_val_adj<1e-300)]=1e-300
myvolcano(F_MCAMvsF_F3,gene.plot = 50,filename = "Figure 7.Fibroblast_diffgene_labeled",
          logfc.cutoff = 1,text.repel = T,width = 26,height = 20,text.size = 1)
myvolcano(F_MCAMvsF_F3,gene.plot = 50,filename = "Figure 7.Fibroblast_diffgene_nolabel",
          logfc.cutoff = 1,text.repel = F,pt.size = 2,width = 5.5)
