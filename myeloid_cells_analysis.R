### myeloid part----
myeloid.subset=subset(immune.combined,idents="Myeloid")
myeloid.subset=subset(Myeloid.subset,subset=nFeature_RNA>400&percent.mt<=20)
colors_myeloid=colors_crc[10:20]

myeloid.list=SplitObject(myeloid.subset,split.by = "patients")
myeloid.list <- lapply(X = myeloid.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = myeloid.list)
length(features)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

library("future.apply")
plan(strategy = "multisession", workers = 4) ## Run in parallel on local computer
options(future.globals.maxSize = 64*1024^3)

myeloid.anchors <- FindIntegrationAnchors(object.list = myeloid.list, anchor.features = features)
myeloid.combined <- IntegrateData(anchorset = myeloid.anchors)
plan(strategy = "multisession", workers = 1)

DefaultAssay(myeloid.combined)="RNA"
myeloid.combined=ScaleData(myeloid.combined)
DefaultAssay(myeloid.combined)="integrated"
myeloid.combined=ScaleData(myeloid.combined)
myeloid.combined=RunPCA(myeloid.combined)

myeloid.combined=FindNeighbors(myeloid.combined,dims = 1:15)
myeloid.combined=FindClusters(myeloid.combined,resolution = 0.3)
myeloid.combined=RunUMAP(myeloid.combined,reduction="pca",dims = 1:15)

## find marker genes ----
DefaultAssay(myeloid.combined)="RNA"
myeloid.markers=myfindmarkers(myeloid.combined,filename = "myeloid",colors = colors_myeloid)
new.cluster.ids <- c("M01_Mono_CD14","M02_Mac_CXCL9","M03_cDC_CD1C","M04_Mac_ZNF331",
                     "M05_Mac_SPP1","M06_Mac_HSPA6","M07_Mac_SERPINB2","M08_cycling_MKI67",
                     "M09_cDC_CPNE3","M10_Mono_FCGR3A","M11_cDC_LAMP3")
names(x = new.cluster.ids) <- levels(x = myeloid.combined)
myeloid.combined <- RenameIdents(object = myeloid.combined, new.cluster.ids)

##figure plot ----
width.ppi=6.5
height.ppi=5
pdf("Figure.Myeloid_dimplot.pdf",width = width.ppi*1.2,height = height.ppi)
DimPlot(object = myeloid.combined, reduction = 'umap',label = ,cols = colors_myeloid,raster = F)
dev.off()

pdf("Figure.Myeloid_dimplot_organs.pdf",width = width.ppi*2.1,height = height.ppi*3)
DimPlot(object = myeloid.combined, reduction = 'umap',label = F,cols = colors_myeloid,raster = F,split.by = "organs",
        ncol = 2,pt.size = 1)
dev.off()

## Figure ----
myeloid.combined=subset(myeloid.combined,subset=patients_organ%in%levels(myeloid.combined$patients_organ)[-c(15,22,29)])
myeloid.combined$patients_organ=paste0(myeloid.combined$organs,"_",myeloid.combined$patients)
myeloid.combined$patients_organ=factor(myeloid.combined$patients_organ,
                                       levels = sort(names(table(myeloid.combined$patients_organ))))
myeloid.combined=subset(myeloid.combined,subset=patients_organ!="CNL_s0920")
percentplot_choose(seurat_object = myeloid.combined,filename = "Figure.myeloid_organ_patients",choose = "patients_organ",
                   width = 11.2,height = 4,legend.key.size = 4,angle = 90,colors = colors_myeloid)

M_CXCL9vsM_SPP1=FindMarkers(myeloid.combined,ident.1 = "M02_Mac_CXCL9",ident.2 = "M05_Mac_SPP1")
M_CXCL9vsM_SPP1$p_val_adj[which(M_CXCL9vsM_SPP1$p_val_adj<1e-300)]=1e-300
myvolcano(M_CXCL9vsM_SPP1,gene.plot = 20,filename = "Figure.Mac_diffgene_labeled",
          logfc.cutoff = 1,text.repel = T,width = 13,height = 10)
myvolcano(M_CXCL9vsM_SPP1,gene.plot = 20,filename = "Figure.Mac_diffgene_nolabel",
          logfc.cutoff = 1,text.repel = F,pt.size = 2)



