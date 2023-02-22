### B part----
load("./immune_delneutro.RData")
B.subset=subset(immune.combined,idents=c("B","plasma"))
B.subset=subset(B.subset,subset=nFeature_RNA>400&percent.mt<=20)
colors_B=colors_crc[1:9]
setwd("./immunecell/Figures/")
dir.create("B")
setwd("./immunecell/Figures/B")

B.list=SplitObject(B.subset,split.by = "patients")
B.list <- lapply(X = B.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = B.list)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

library("future.apply")
plan(strategy = "multisession", workers = 4) ## Run in parallel on local computer
options(future.globals.maxSize = 64*1024^3)
B.anchors <- FindIntegrationAnchors(object.list = B.list, anchor.features = features)
B.combined <- IntegrateData(anchorset = B.anchors)
plan(strategy = "multisession", workers = 1)

DefaultAssay(B.combined)="RNA"
B.combined=ScaleData(B.combined)
DefaultAssay(B.combined)="integrated"
B.combined=ScaleData(B.combined)

B.combined=RunPCA(B.combined)
ElbowPlot(B.combined,ndims = 50)
B.combined=FindNeighbors(B.combined,dims = 1:15)
B.combined=FindClusters(B.combined,resolution = 0.3)
B.combined=RunUMAP(B.combined,reduction="pca",dims = 1:15)

## find marker genes ----
setwd("./immunecell/Figures/B")
DefaultAssay(B.combined)="RNA"
DimPlot(B.combined,label=T,cols = colors_B)
B.markers=myfindmarkers(B.combined,filename = "B",colors = colors_B)
new.cluster.ids <- c("B01_plasma_IGKC","B02_plasma_IGLC3","B03_B_MS4A1","B01_plasma_IGKC",
                     "B04_plasma_IGHG2","B05_plasma_IGLL1","B06_B_IGHD","B07_plasma_PPIB",
                     "B08_plasma_IGLJ3","B09_GCB_RGS13")
names(x = new.cluster.ids) <- levels(x = B.combined)
B.combined <- RenameIdents(object = B.combined, new.cluster.ids)
