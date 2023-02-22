##tumorcell part----
##subset tumorcell----
colors_tumorcell=colors_tumor[1:11]
tumorcell.subset=subset(tumor.combined,idents="Tumorcell")

tumorcell.list=SplitObject(tumorcell.subset,split.by = "patients")
tumorcell.list <- lapply(X = tumorcell.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = tumorcell.list)
length(features)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

library("future.apply")
plan(strategy = "multisession", workers = 4) ## Run in parallel on local computer
options(future.globals.maxSize = 64*1024^3)
tumorcell.anchors <- FindIntegrationAnchors(object.list = tumorcell.list, anchor.features = features)
tumorcell.combined <- IntegrateData(anchorset = tumorcell.anchors)
plan(strategy = "multisession", workers = 1)
DefaultAssay(tumorcell.combined)="RNA"
tumorcell.combined=ScaleData(tumorcell.combined)
DefaultAssay(tumorcell.combined)="integrated"
tumorcell.combined=ScaleData(tumorcell.combined)
tumorcell.combined=RunPCA(tumorcell.combined)
ElbowPlot(tumorcell.combined,ndims = 50)
tumorcell.combined=FindNeighbors(tumorcell.combined,dims = 1:20)
tumorcell.combined=FindClusters(tumorcell.combined,resolution = 0.4)
tumorcell.combined=RunUMAP(tumorcell.combined,reduction="pca",dims = 1:20)

## find marker genes ----
DefaultAssay(tumorcell.combined)="RNA"
tumorcell.combined=subset(tumorcell.combined,subset=percent.mt<=40)
tumorcell.markers=myfindmarkers(tumorcell.combined,filename = "tumorcell",colors = colors_tumorcell)

new.cluster.ids <- c("Tu01_AREG","Tu02_DEFA5","Tu03_SRRM2","Tu04_RGMB",
                     "Tu05_PCNA","Tu06_NKD1","Tu07_MKI67","Tu08_GNG13",
                     "Tu09_MUC2","Tu10_COL3A1","Tu11_PLA2G2A","Tu01_AREG")
names(x = new.cluster.ids) <- levels(x = tumorcell.combined)
tumorcell.combined <- RenameIdents(object = tumorcell.combined, new.cluster.ids)

## Figure plot ----
width.ppi=6.5
height.ppi=5
setwd("./nonimmune/tumorcell")

## Figure ----
pdf("Figure.tumorcell_dimplot.pdf",width = width.ppi*1.1,height = height.ppi)
DimPlot(object = tumorcell.combined, reduction = 'umap',label = F,cols = colors_tumorcell)
dev.off()

pdf("Figure.tumorcell_organs_dimplot.pdf",width = width.ppi*1.8,height = height.ppi)
DimPlot(object = tumorcell.combined, reduction = 'umap',label = F,cols = colors_tumorcell,raster = F,split.by = "organs",
        ncol = 2)
dev.off()

tumorcell.markers=myfindmarkers(tumorcell.combined,filename = "Figure.tumorcell",colors = colors_tumorcell)
## Figure ----
stacked_violin_plot(gene = c("LGR5","EPCAM","CDH1"),col = colors_tumorcell,filename = "Figure.stem_gene",
                    seurat_object = tumorcell.combined,Mean = F,flip = F,width = 8,height = 6)
