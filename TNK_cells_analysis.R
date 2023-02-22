##TNK part----
tnk.subset=subset(immune.combined,idents="TNK")
tnk.subset=subset(tnk.subset,subset=percent.mt<=20&nFeature_RNA>200)

DefaultAssay(tnk.subset)="integrated"
tnk.subset=ScaleData(tnk.subset)
tnk.subset=RunPCA(tnk.subset)
ElbowPlot(tnk.subset,ndims = 50)
tnk.subset=FindNeighbors(tnk.subset,dims = 1:15)
tnk.subset=FindClusters(tnk.subset,resolution = 1)
tnk.subset=RunUMAP(tnk.subset,reduction="pca",dims = 1:15)

DefaultAssay(tnk.subset)="RNA"
pdf("tnk_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = tnk.subset, reduction = 'umap',label = T,cols = colors_crc)
dev.off()

tnk.subset$patients_organ=paste0(tnk.subset$organs,"_",tnk.subset$patients)
tnk.subset$patients_organ=factor(tnk.subset$patients_organ,
                                 levels = sort(names(table(tnk.subset$patients_organ))))
##delete only 5 cells sample----
tnk.subset=subset(tnk.subset,subset=patients_organ!="CNL_s0920")
#rm(immune.combined)
#gc()
antibody=c("CD3.pAbO","CD4.pAbO","CD8.pAbO","CD279.pAbO","TIM3.pAbO","LAG3.pAbO","CD25.pAbO",
           "GITR.pAbO","CD69.pAbO","CD103.pAbO","CCR7.pAbO","CD127.pAbO","CXCR3.pAbO","CD7.pAbO","CD27.pAbO","CD56.pAbO","CD335.pAbO",
           "TCRgd.pAbO","CD14.pAbO","CD16.pAbO","CD11b.pAbO","CD11c.pAbO","CD273.pAbO","CD274.pAbO",
           "CD19.pAbO","CD20.pAbO","CD24.pAbO","CD38.pAbO","CXCR5.pAbO","IgD.pAbO","IgG.pAbO",
           "CD45.pAbO","CD45RA.pAbO","CD45RO.pAbO")
dir.create("TNK")
setwd("/GPUFS/scut_zxlian_4/longjie/20220111_CRC/20220205_results/immunecell/Figures/TNK")
DefaultAssay(tnk.subset)="antibody"
myfeatureplot(tnk.subset,filename = "tnk_antibody",genename = antibody,reduction = "umap")

new.cluster.ids <- c("CD4","CD8","CD8","CD4","CD8","CD4","NK01_NCAM1","CD8","CD8","CD8","NK02_FCGR3A",
                     "CD4","gdT","CD4","cycling_TNK","CD4","CD4","CD4","CD8","CD4",
                     "NK02_FCGR3A","NK01_NCAM1")
names(x = new.cluster.ids) <- levels(x = tnk.subset)
tnk.subset <- RenameIdents(object = tnk.subset, new.cluster.ids)
save(tnk.subset,file = "TNK.RData")

##cd8 part----
##subset cd8----
load("TNK.RData")
colors_cd8=colors_crc[27:34]
DefaultAssay(tnk.subset)="RNA"
cd8.subset=subset(tnk.subset,idents="CD8")

cd8.list=SplitObject(cd8.subset,split.by = "patients")
cd8.list <- lapply(X = cd8.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = cd8.list)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

library("future.apply")
plan(strategy = "multisession", workers = 4) ## Run in parallel on local computer
options(future.globals.maxSize = 16*1024^3)

cd8.anchors <- FindIntegrationAnchors(object.list = cd8.list, anchor.features = features)
cd8.combined <- IntegrateData(anchorset = cd8.anchors)
DefaultAssay(cd8.combined)="RNA"
cd8.combined=ScaleData(cd8.combined)
DefaultAssay(cd8.combined)="integrated"
cd8.combined=ScaleData(cd8.combined)
cd8.combined=RunPCA(cd8.combined)
cd8.combined=FindNeighbors(cd8.combined,dims = 1:20)
cd8.combined=FindClusters(cd8.combined,resolution = 0.4)
cd8.combined=RunUMAP(cd8.combined,reduction="pca",dims = 1:20)

## find marker genes ----
DefaultAssay(cd8.combined)="RNA"
cd8.markers=myfindmarkers(cd8.combined,filename = "CD8",colors = colors)
new.cluster.ids <- c("T01_CD8_RGCC","T02_CD8_CCL4","T03_CD8_CXCL13","T04_CD8_GZMK","T05_CD8_HSPA6",
                     "T06_CD8_CX3CR1","T07_CD8_ITGA1")
names(x = new.cluster.ids) <- levels(x = cd8.combined)
cd8.combined <- RenameIdents(object = cd8.combined, new.cluster.ids)

pdf("Figure.cd8_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = cd8.combined, reduction = 'umap',label = F,cols = colors_cd8)
dev.off()

pdf("Figure.cd8_dimplot_samples.pdf",width = width.ppi*1.1,height = height.ppi*1.5)
DimPlot(object = cd8.combined, reduction = 'umap',label = F,cols = colors_cd8,split.by = "organs",
        ncol = 2)
dev.off()

##CD4 part----
##subset cd4----
cd4.subset=subset(tnk.subset,idents="CD4")

cd4.list=SplitObject(cd4.subset,split.by = "patients")
cd4.list <- lapply(X = cd4.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = cd4.list)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

library("future.apply")
plan(strategy = "multisession", workers = 4) ## Run in parallel on local computer
options(future.globals.maxSize = 64*1024^3)

cd4.anchors <- FindIntegrationAnchors(object.list = cd4.list, anchor.features = features)
cd4.combined <- IntegrateData(anchorset = cd4.anchors)
plan(strategy = "multisession", workers = 1)
DefaultAssay(cd4.combined)="RNA"
cd4.combined=ScaleData(cd4.combined)
DefaultAssay(cd4.combined)="integrated"
cd4.combined=ScaleData(cd4.combined)
cd4.combined=RunPCA(cd4.combined)
cd4.combined=FindNeighbors(cd4.combined,dims = 1:12)
cd4.combined=FindClusters(cd4.combined,resolution = 0.45)
cd4.combined=RunUMAP(cd4.combined,reduction="pca",dims = 1:12)

DefaultAssay(cd4.combined)="RNA"
cd4.markers=myfindmarkers(cd4.combined,filename = "CD4",colors = colors_cd4)

new.cluster.ids <- c("T09_CD4_BHLHE40","T10_CD4_SELL","T11_CD4_HSPA6","T12_Treg_FOXP3",
                     "T13_CD4_ZBTB10","T14_Treg_TNFRSF9","T15_CD4_CXCL13")
names(x = new.cluster.ids) <- levels(x = cd4.combined)
cd4.combined <- RenameIdents(object = cd4.combined, new.cluster.ids)

width.ppi=6.5
height.ppi=5
pdf("Figure.CD4_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = cd4.combined, reduction = 'umap',label = F,cols = colors_cd4)
dev.off()

pdf("Figure.CD4_dimplot_samples.pdf",width = width.ppi*1.1,height = height.ppi*1.5)
DimPlot(object = cd4.combined, reduction = 'umap',label = F,cols = colors_cd4,split.by = "organs",
        ncol = 2)
dev.off()
