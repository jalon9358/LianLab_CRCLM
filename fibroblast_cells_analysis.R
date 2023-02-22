##fibroblast part----
##subset fibroblast----
fibroblast.subset=subset(tumor.combined,idents="Fibroblast")
colors_fibroblast=colors_tumor[12:17]

fibroblast.list=SplitObject(fibroblast.subset,split.by = "patients")
fibroblast.list <- lapply(X = fibroblast.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = fibroblast.list)
length(features)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

library("future.apply")
plan(strategy = "multisession", workers = 4) ## Run in parallel on local computer
options(future.globals.maxSize = 64*1024^3)
fibroblast.anchors <- FindIntegrationAnchors(object.list = fibroblast.list, anchor.features = features)
fibroblast.combined <- IntegrateData(anchorset = fibroblast.anchors)
plan(strategy = "multisession", workers = 1)

DefaultAssay(fibroblast.combined)="RNA"
fibroblast.combined=ScaleData(fibroblast.combined)
DefaultAssay(fibroblast.combined)="integrated"
fibroblast.combined=ScaleData(fibroblast.combined)
fibroblast.combined=RunPCA(fibroblast.combined)
ElbowPlot(fibroblast.combined,ndims = 50)
fibroblast.combined=FindNeighbors(fibroblast.combined,dims = 1:12)
fibroblast.combined=FindClusters(fibroblast.combined,resolution = 0.6)
fibroblast.combined=RunUMAP(fibroblast.combined,reduction="pca",dims = 1:12)

## find marker genes ----
DefaultAssay(fibroblast.combined)="RNA"
fibroblast.markers=myfindmarkers(fibroblast.combined,filename = "fibroblast",colors = colors_fibroblast)
new.cluster.ids <- c("F01_fibroblast_PRELP","F02_fibrblast_MCAM","F03_fibroblast_CXCL14",
                     "F04_fibroblast_C3","F05_fibroblast_COCH","F06_cycling_MKI67")
names(x = new.cluster.ids) <- levels(x = fibroblast.combined)
fibroblast.combined <- RenameIdents(object = fibroblast.combined, new.cluster.ids)
fibroblast.combined$cluster=fibroblast.combined@active.ident

## delete mitochondial gene high cells ----
fibroblast.combined=subset(fibroblast.combined,subset=percent.mt<30)

## Figure ----
width.ppi=6.5
height.ppi=5
pdf("Figure.fibroblast_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = fibroblast.combined, reduction = 'umap',label = F,cols = colors_fibroblast)
dev.off()

pdf("Figure.fibroblast_dimplot_organs.pdf",width = width.ppi*1.8,height = height.ppi*1)
DimPlot(object = fibroblast.combined, reduction = 'umap',label = F,cols=colors_fibroblast,
        split.by = "organs",
        ncol = 2)
dev.off()

## Figure ----
DefaultAssay(fibroblast.combined)="RNA"
fibroblast.markers=myfindmarkers(fibroblast.combined,filename = "Figure.fibroblast",colors = colors_fibroblast)
