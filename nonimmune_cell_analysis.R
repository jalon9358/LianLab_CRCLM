setwd("./nonimmune")
load("tumor_integrated.RData")

DefaultAssay(tumor.combined)="RNA"
tumor.combined=ScaleData(tumor.combined)
DefaultAssay(tumor.combined)="integrated"
tumor.combined=ScaleData(tumor.combined)
tumor.combined=RunPCA(tumor.combined)
ElbowPlot(tumor.combined,ndims = 50)
tumor.combined=FindNeighbors(tumor.combined,dims = 1:20)
tumor.combined=FindClusters(tumor.combined,resolution = 0.1)
tumor.combined=RunUMAP(tumor.combined,reduction="pca",dims = 1:20)

##scrublet delete doublet----
##***********************###
tumor.combined$batch=tumor.combined$patients
tumor.batch=SplitObject(tumor.combined,split.by = "batch")
dir.create("scrublet")
setwd("./nonimmune/scrublet")
for (i in 1:length(tumor.batch)) {
  count.matrix=tumor.batch[[i]]@assays$RNA@counts
  count.matrix=t(count.matrix)
  write.csv(count.matrix,file = paste0(names(tumor.batch)[i],".csv"))
  rm(count.matrix)
  gc()
}

## python scrublet.py----
directory=dir("./nonimmune/doublet_score")
path="./nonimmune/doublet_score"

for (i in directory) {
  score_split=read.csv(paste0(path,"/",i),header = T,row.names = 1)
  if(i==directory[1]){
    score=score_split
  }else{
    score=rbind(score,score_split)
  }
}

score=score[colnames(tumor.combined),]
tumor.combined$doublet.score=score$score
tumor.combined$predicted.doublet=score$prediction
table(tumor.combined$predicted.doublet)
tumor.combined$doublet=ifelse(tumor.combined$predicted.doublet=="True","doublet","singlet")
tumor.combined$cluster=tumor.combined@active.ident
Idents(tumor.combined)="cluster"
myfeatureplot(tumor.combined,genename = "doublet.score",raster = F,
              filename = "doublet",reduction = "umap",pt.size = 0.5)

Idents(tumor.combined)="doublet"
p1=VlnPlot(tumor.combined,features = "nFeature_RNA",pt.size = 0)
p2=VlnPlot(tumor.combined,features = "nCount_RNA",pt.size = 0)
pdf("doublet_geneno.pdf",width = width.ppi*1.5,height = height.ppi)
multiplot(p1,p2,cols = 2)
dev.off()

setwd("./nonimmune/")
DefaultAssay(tumor.combined)="RNA"
tumor.markers=myfindmarkers(tumor.combined,filename = "tumor",colors = colors2)

new.cluster.ids <- c("Tumorcell","Tumorcell","Fibroblast","Fibroblast","Tumorcell","Endo",
                     "Tumorcell","Tumorcell","Tumorcell")  
names(x = new.cluster.ids) <- levels(x = tumor.combined)
tumor.combined <- RenameIdents(object = tumor.combined, new.cluster.ids)
