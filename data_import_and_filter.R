#------Import BD data------
#1.1.source functions and libraries----
source("/Users/jalon/Desktop/newest scripts/functions/source.R")
source("/Users/jalon/Desktop/newest scripts/functions/read_BDdata.R")
#1.2.import BD datasheet----
path="/Users/jalon/Desktop/analysis/wangfei_CRC_sc/rawdata"
merge_data=read_BDdata(path = path)

for (i in 1:length(merge_data)) {
  rownames(merge_data[[i]])=paste0(merge_data[[i]][,"sampleinfo"],"-",
                                   str_split_fixed(rownames(merge_data[[i]]),"-",2)[,2])
}

for (i in 1:length(merge_data)) {
  cat("sample ",i,"\n")
  merge_data[[i]]=subset(merge_data[[i]],select=-c(sample_name,sample_tag,sampleinfo))
  merge_data[[i]]=as.sparse(merge_data[[i]])
}

for (i in 1:length(merge_data)) {
  merge_data[[i]]=t(merge_data[[i]])
}

for (i in 1:length(merge_data)) {
  colnames(merge_data[[i]])=paste0(colnames(merge_data[[i]]),"-",str_split_fixed(names(merge_data)[i],"_",2)[,2])
}

##1.3.get antibody name----
antibody_name=c()
for(i in 1:length(merge_data)){
  rownames(merge_data[[i]])[which(grepl(pattern = "\\.pAbO$",
                                        rownames(merge_data[[i]])))]=paste0(str_split_fixed(grep(pattern = "\\.pAbO$",
                                                                                                 rownames(merge_data[[i]]),value = T),"\\.",2)[,1],".pAbO")
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="TCR.pAbO")]="TCRgd.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="LAG.pAbO")]="LAG3.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="Lag3.pAbO")]="LAG3.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="Tim3.pAbO")]="TIM3.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="CD223.pAbO")]="Lag3.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="CD357.pAbO")]="GITR.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="CD185.pAbO")]="CXCR5.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="CD183.pAbO")]="CXCR3.pAbO"
  rownames(merge_data[[i]])[which(rownames(merge_data[[i]])=="CD197.pAbO")]="CCR7.pAbO"
  if(i==1){
    antibody_name=grep(pattern = "\\.pAbO$",rownames(merge_data[[i]]),value = T)
  }else{
    antibody_name=union(antibody_name,grep(pattern = "\\.pAbO$",rownames(merge_data[[i]]),value = T))
  }
}

##1.4.create immune and tumor object----
immune_merge=c(list(merge_data$s0920_3_4),list(merge_data$s0920_5_7),
               list(merge_data$s1125_3_4),list(merge_data$s1125_1_2),
               list(merge_data$s1231_3_4),list(merge_data$s1231_5_6),
               list(merge_data$s0107_3_4),list(merge_data$s0107_5_7),
               list(merge_data$s0115_3_4),list(merge_data$s0115_5_7),
               list(merge_data$s0813_3_4),list(merge_data$s0813_5_7),
               list(merge_data$s0816_2_4))

tumor_merge=c(list(merge_data$s1231_1_2),list(merge_data$s0920_1_2),list(merge_data$s0813_1_2),
              list(merge_data$s0115_1_2),list(merge_data$s0107_1_2),list(merge_data$s0816_1))


###immunecell part----
gene_inter=c()
gene_union=c()
for(i in 1:length(immune_merge)){
  gene=rownames(immune_merge[[i]])
  if(i==1){
    gene_inter=gene
    gene_union=gene
  }else{
    gene_inter=intersect(gene_inter,gene)
    gene_union=union(gene_union,gene)
  }
}

##extend the matrix according to union genes of all samples----
for (i in 1:length(immune_merge)) {
  gene_insert=setdiff(gene_union,rownames(immune_merge[[i]]))
  insert_matrix=as.data.frame(matrix(0,nrow = length(gene_insert),ncol = ncol(immune_merge[[i]])))
  row.names(insert_matrix)=gene_insert
  colnames(insert_matrix)=colnames(immune_merge[[i]])
  norow=nrow(immune_merge[[i]])
  insert_matrix=as.sparse(insert_matrix)
  immune_merge[[i]]=rbind(immune_merge[[i]],insert_matrix)
  immune_merge[[i]]=immune_merge[[i]][sort(gene_union),]
  cat("sample ",i, " insert gene number ",nrow(insert_matrix),"\n")
}

#merge all gene filtered data----
for (i in 1:length(immune_merge)) {
  cat("sample ",i,"\n")
  if(i==1){
    immune_count=immune_merge[[i]]
  }else{
    immune_count=cbind(immune_count,immune_merge[[i]])
  }
}

noncodegene=grep(pattern = "\\.AS|^LINC|^MIR[0-9]*|\\.[I|O]T1|\\.DT$|^[A-Z][A-Z][0-9]{3,10}\\.[0-9]{1,3}|^[A-Z][0-9]{3,10}\\.[0-9]{1,3}",gene_union,value = T)
immune_count=immune_count[setdiff(rownames(immune_count),noncodegene),]

#split RNA and antibody data----
antibody_matrix=immune_count[grep(pattern = "\\.pAbO$",rownames(immune_count)),]
RNA_matrix=immune_count[-grep(pattern = "\\.pAbO$",rownames(immune_count)),]

###create immune seurat object----
immune_object=CreateSeuratObject(RNA_matrix,min.cells = 50)
immune_object[["antibody"]] <- CreateAssayObject(counts = antibody_matrix)

###2.add samples information----
immune_object$patients=immune_object$orig.ident
sampleinfo_split=str_split_fixed(colnames(immune_object),"-",3)
immune_object$sampletag=paste0(sampleinfo_split[,1],"-",sampleinfo_split[,3])
sampleinfo=read.csv("/Users/jalon/Desktop/analysis/wangfei_CRC_sc/sampleinfo.csv")
sampleinfo$sampleinfo=paste0(sampleinfo$date,"_",sampleinfo$sampletag,"-",sampleinfo$no)
LCT_samples=subset(sampleinfo,subset=cell_type=="LCT")[,"sampleinfo"]
CCT_samples=subset(sampleinfo,subset=cell_type=="CCT")[,"sampleinfo"]
CCL_samples=subset(sampleinfo,subset=cell_type=="CCL")[,"sampleinfo"]
CNL_samples=subset(sampleinfo,subset=cell_type=="CNL")[,"sampleinfo"]
LCL_samples=subset(sampleinfo,subset=cell_type=="LCL")[,"sampleinfo"]
LNL_samples=subset(sampleinfo,subset=cell_type=="LNL")[,"sampleinfo"]
PBL_samples=subset(sampleinfo,subset=cell_type=="PBL")[,"sampleinfo"]
immune_object$organs=ifelse(immune_object$sampletag%in%CCL_samples,"CCL",
                            ifelse(immune_object$sampletag%in%CNL_samples,"CNL",
                                   ifelse(immune_object$sampletag%in%LCL_samples,"LCL",
                                          ifelse(immune_object$sampletag%in%LNL_samples,"LNL",
                                                 ifelse(immune_object$sampletag%in%PBL_samples,"PBL",
                                                        ifelse(immune_object$sampletag%in%LCT_samples,"LCT",
                                                               ifelse(immune_object$sampletag%in%CCT_samples,"CCT","Others")))))))

### 3.cell filter----
immune_object[["percent.mt"]] <- PercentageFeatureSet(immune_object, pattern = "^MT\\.")
immune_object[["percent.ribo"]] <- PercentageFeatureSet(immune_object, pattern = "^RP[S|L]")
#calculate log10GenesPerUMI
immune_object$log10GenesPerUMI <- log10(immune_object$nFeature_RNA)/log10(immune_object$nCount_RNA)

width.ppi=6.5
height.ppi=5
pdf("cutoff.pdf", width = width.ppi*1, height = height.ppi*1.5)
p1=VlnPlot(immune_object,features = c("nFeature_RNA"),pt.size = 0)+geom_hline(yintercept = 200)+
  geom_hline(yintercept = 4000)+NoLegend()+
  labs(subtitle=paste0("depleted cellno is ",length(which(immune_object$nFeature_RNA>4000|immune_object$nFeature_RNA<200))))
p2=VlnPlot(immune_object,features = c("nCount_RNA"),pt.size = 0)+geom_hline(yintercept = 15000)+NoLegend()+
  labs(subtitle=paste0("depleted cellno is ",length(which(immune_object$nCount_RNA>15000))))
p3=VlnPlot(immune_object,features = c("percent.mt"),pt.size = 0)+geom_hline(yintercept = 25)+NoLegend()+
  labs(subtitle=paste0("depleted cellno is ",length(which(immune_object$percent.mt>25))))
multiplot(plotlist = list(p1,p2,p3),cols = 2)
dev.off()
immune_object=subset(immune_object,subset=nFeature_RNA>200 & nFeature_RNA<4000)
immune_object=subset(immune_object,subset=percent.mt<25)
immune_object=subset(immune_object,subset=nCount_RNA<15000)
colname_split=str_split_fixed(colnames(immune_object),"-",3)
immune_object$batch=paste0(str_split_fixed(colname_split[,1],"_",2)[,1],"_",colname_split[,3])

### 4.batch effect correct----
#bactch effect correct----
immune_object$patients=str_split_fixed(colnames(immune_object),"_",2)[,1]
immune.list=SplitObject(immune_object,split.by = "patients")

immune.list <- lapply(X = immune.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
})

features <- SelectIntegrationFeatures(object.list = immune.list,nfeatures = 2500)
if(length(grep(pattern = "^RP[L|S]|^MT\\.",features))!=0){
  features=features[-grep(pattern = "^RP[L|S]|^MT\\.",features)]
}

immune.anchors <- FindIntegrationAnchors(object.list = immune.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

### 5.dimension reduction and clustering----
DefaultAssay(immune.combined)="RNA"
immune.combined=ScaleData(immune.combined)
DefaultAssay(immune.combined)="integrated"
immune.combined=ScaleData(immune.combined)
immune.combined=RunPCA(immune.combined)
immune.combined=FindNeighbors(immune.combined,dims = 1:20)
immune.combined=FindClusters(immune.combined,resolution = 0.5)
immune.combined=RunUMAP(immune.combined,reduction="pca",dims = 1:20)
set.seed(1314)
colors_crc=sample(colors_aa,49)

pdf("total_immunecell.pdf",width = width.ppi*1.3,height = height.ppi*1.2)
DimPlot(object = immune.combined, reduction = 'umap',label = T,cols = colors_crc,raster = F)
dev.off()

pdf("total_immunecell_organs.pdf",width = width.ppi*2.1,height = height.ppi*3)
DimPlot(object = immune.combined, reduction = 'umap',label = T,cols = colors_crc,raster = F,
        split.by = "organs",ncol = 2)
dev.off()

topgene=myfindmarkers(immune.combined,filename = "immune_total",colors = colors_crc)

antibody=rownames(immune.combined@assays$antibody)
antibody=c("CD3.pAbO","CD4.pAbO","CD8.pAbO","CD279.pAbO","TIM3.pAbO","LAG3.pAbO","CD25.pAbO",
           "GITR.pAbO","CD69.pAbO","CD103.pAbO","CCR7.pAbO","CD127.pAbO","CXCR3.pAbO","CD7.pAbO","CD27.pAbO","CD56.pAbO","CD335.pAbO",
           "TCRgd.pAbO","CD14.pAbO","CD16.pAbO","CD11b.pAbO","CD11c.pAbO","CD273.pAbO","CD274.pAbO",
           "CD19.pAbO","CD20.pAbO","CD24.pAbO","CD38.pAbO","CXCR5.pAbO","IgD.pAbO","IgG.pAbO",
           "CD45.pAbO","CD45RA.pAbO","CD45RO.pAbO")
gene_antibody=c("CD3D","CD4","CD8A","PDCD1","HAVCR2","LAG3","IL2RA","TNFRSF18","CD69",
                "ITGAE","CCR7","IL7R","CXCR3","CD7","CD27","NCAM1","NCR1","TRGC1","CD14","FCGR3A",
                "ITGAM","ITGAX","PDCD1LG2","CD274","CD19","MS4A1","CD24","CD38","CXCR5","IGHD",
                "IGHG1","IGHG2","IGHG3","IGHG4","PTPRC","PTPRC","PTPRC")

DefaultAssay(immune.combined)="antibody"
immune.combined=NormalizeData(immune.combined,assay = "antibody",normalization.method = "CLR")
myfeatureplot(immune.combined,filename = "antibody",genename = antibody,reduction = "umap")

DefaultAssay(immune.combined)="RNA"
myfeatureplot(immune.combined,filename = "RNA",genename = gene_antibody,reduction = "umap")
immune.combined$samples=immune.combined$organs
percentplot_sample(immune.combined,filename = "Immune",colors = colors_crc)

DefaultAssay(immune.combined)="antibody"
immune.combined=NormalizeData(immune.combined,assay = "antibody",normalization.method = "CLR")

### 6.scrublet delete doublet ----
immune.batch=SplitObject(immune.combined,split.by = "batch")
for (i in 1:length(immune.batch)) {
  count.matrix=immune.batch[[i]]@assays$RNA@counts
  count.matrix=t(count.matrix)
  write.csv(count.matrix,file = paste0(names(immune.batch)[i],".csv"))
  rm(count.matrix)
}

directory=dir("./immunecell/doublet_score")
path="./immunecell/doublet_score"

for (i in directory) {
  score_split=read.csv(paste0(path,"/",i),header = T,row.names = 1)
  if(i==directory[1]){
    score=score_split
  }else{
    score=rbind(score,score_split)
  }
}
score=score[colnames(immune.combined),]
immune.combined$doublet.score=score$score
immune.combined$predicted.doublet=score$prediction
immune.combined$cluster=immune.combined@active.ident
immune.combined$doublet=ifelse(immune.combined$predicted.doublet=="True","doublet","singlet")
DimPlot(immune.combined,split.by = "doublet")

Idents(immune.combined)="cluster"
myfeatureplot(immune.combined,genename = "doublet.score",
              filename = "doublet",reduction = "umap",pt.size = 0.5)

Idents(immune.combined)="doublet"
p1=VlnPlot(immune.combined,features = "nFeature_RNA",pt.size = 0)
p2=VlnPlot(immune.combined,features = "nCount_RNA",pt.size = 0)
width.ppi=6.5
height.ppi=5
pdf("doublet_geneno.pdf",width = width.ppi*1.5,height = height.ppi)
multiplot(p1,p2,cols = 2)
dev.off()
immune.combined=subset(immune.combined,subset=doublet=="singlet")
Idents(immune.combined)="cluster"

## 7. Figures plot ----
##delete neutrophils----
immune.combined=subset(immune.combined,idents=levels(immune.combined)[-14])
DefaultAssay(immune.combined)="RNA"
immune.combined=ScaleData(immune.combined)
DefaultAssay(immune.combined)="integrated"
immune.combined=ScaleData(immune.combined)
immune.combined=RunPCA(immune.combined)
ElbowPlot(immune.combined,ndims = 50)
immune.combined=FindNeighbors(immune.combined,dims = 1:20)
immune.combined=FindClusters(immune.combined,resolution = 0.5)
immune.combined=RunUMAP(immune.combined,reduction="pca",dims = 1:20)

width.ppi=6.5
height.ppi=5

pdf("total_immunecell.pdf",width = width.ppi*1.3,height = height.ppi*1.2)
DimPlot(object = immune.combined, reduction = 'umap',label = T,cols = colors_crc,raster = F)
dev.off()

pdf("total_immunecell_organs.pdf",width = width.ppi*2.1,height = height.ppi*3)
DimPlot(object = immune.combined, reduction = 'umap',label = T,cols = colors_crc,raster = F,
        split.by = "organs",ncol = 2)
dev.off()

pdf("total_immunecell_patients.pdf",width = width.ppi*3.1,height = height.ppi*3)
DimPlot(object = immune.combined, reduction = 'umap',label = T,cols = colors_crc,raster = F,
        split.by = "patients",ncol = 3)
dev.off()

topgene=myfindmarkers(immune.combined,filename = "immune_total",colors = colors_crc)

new.cluster.ids <- c("TNK","TNK","plasma","TNK","plasma","TNK","TNK","Myeloid","TNK","B","TNK","TNK","TNK",
                     "plasma","Myeloid","Myeloid","TNK","Mast cell","Myeloid","plasma","pDC","TNK","Myeloid",
                     "Tumor cells")
names(x = new.cluster.ids) <- levels(x = immune.combined)
immune.combined <- RenameIdents(object = immune.combined, new.cluster.ids)













