##endothelial part----
##subset endothelial----
colors_endothelial=colors_tumor[18:23]
endothelial.subset=subset(tumor.combined,idents="Endo")

endothelial.list=SplitObject(endothelial.subset,split.by = "patients")
endothelial.list <- lapply(X = endothelial.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = endothelial.list)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

library("future.apply")
plan(strategy = "multisession", workers = 4) ## Run in parallel on local computer
options(future.globals.maxSize = 64*1024^3)

endothelial.anchors <- FindIntegrationAnchors(object.list = endothelial.list, anchor.features = features)
endothelial.combined <- IntegrateData(anchorset = endothelial.anchors)
plan(strategy = "multisession", workers = 1)
DefaultAssay(endothelial.combined)="RNA"
endothelial.combined=ScaleData(endothelial.combined)
DefaultAssay(endothelial.combined)="integrated"
endothelial.combined=ScaleData(endothelial.combined)
endothelial.combined=RunPCA(endothelial.combined)
ElbowPlot(endothelial.combined,ndims = 50)
endothelial.combined=FindNeighbors(endothelial.combined,dims = 1:10)
endothelial.combined=FindClusters(endothelial.combined,resolution = 0.2)
endothelial.combined=RunUMAP(endothelial.combined,reduction="pca",dims = 1:10)

## find marker genes ----
DefaultAssay(endothelial.combined)="RNA"
endothelial.markers=myfindmarkers(endothelial.combined,filename = "endothelial",colors = colors_endothelial)

## delete high percent.mt cells ----
endothelial.combined=subset(endothelial.combined,subset=percent.mt<40)
DefaultAssay(endothelial.combined)="RNA"
endothelial.markers=myfindmarkers(endothelial.combined,filename = "endothelial",colors = colors_endothelial)
new.cluster.ids <- c("E01_endothelial_SELP","E02_endothelial_DLL4","E03_endothelial_NOTCH3",
                     "E04_endothelial_CD36","E05_cycling_MKI67","E06_endothelial_CLEC4G")
names(x = new.cluster.ids) <- levels(x = endothelial.combined)
endothelial.combined <- RenameIdents(object = endothelial.combined, new.cluster.ids)
