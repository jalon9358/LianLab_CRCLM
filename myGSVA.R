#****Long Jie,2021-03-01****#
#This function is used to perform GSVA analysis
#@seurat_object:a seurat_object
#@gspath:the correct path of your geneset file
if(!require(GSVA)){BiocManager::install("GSVA");library(GSVA)}
if(!require(GSEABase)){BiocManager::install("GSEABase");library(GSEABase)}
if(!require(pheatmap)){install.packages("pheatmap");library(pheatmap)}
if(!require(limma)){install.packages("limma");library(limma)}
if(!require(biomaRt)){BiocManager::install("biomaRt");library(biomaRt)}

create_supercell_matrix=function(seurat_object,ncells=50){
  dm=as.data.frame(seurat_object@assays$RNA@data)
  split.clusters=split(as.data.frame(t(dm)),seurat_object@active.ident)
  supercell_matrix=data.frame()
  for (i in 1:length(split.clusters)) {
    cat("cluster ",i,"\n")
    split.matrix=split.clusters[[i]]
    along=1:nrow(split.matrix)
    supercells=aggregate(split.matrix,list(ceiling(along/ncells)),mean)
    supercells$clusters=rep(names(split.clusters)[[i]],nrow(supercells))
    supercell_matrix=rbind(supercell_matrix,supercells)
  }
  cluster=supercell_matrix$clusters
  supercell_matrix=subset(supercell_matrix,select=-c(Group.1,clusters))
  supercell_matrix=as.data.frame(t(supercell_matrix))
  return(list(supercell_matrix=supercell_matrix,cluster=cluster))
}

GSVA_cluster_diff=function(seurat_object=NULL,method="ssgsea",
                           gspath="/mnt/ssd1/longjie/newest_scripts/database/gsea/h.all.v7.1.symbols.gmt",
                           species="human",min.sz=2,
                           parallel.sz=12,supercell_no=50,adjPvalueCutoff=0.05,filename="",
                           logFCcutoff=0.1,text.size=2,topn=20,clusters=NULL,pathname=NULL,
                           width = 6,height = 5,
                           cluster_rows=F,cluster_cols=F,fontsize_row=6,
                           col=rev(RColorBrewer::brewer.pal(n=11,name = "RdBu"))){
  
  if(is.null(seurat_object)){
    stop("No seurat object supplied!")
  }else{
    library(Seurat)
    #seurat_object=FindVariableFeatures(seurat_object)
    #seurat_object=subset(seurat_object,features=seurat_object@assays$RNA@var.features)
    print("step1. Create supercells...")
    total_matrix=create_supercell_matrix(seurat_object = seurat_object,ncells = supercell_no)
    supercell_cluster=total_matrix$cluster
    total_matrix=total_matrix$supercell_matrix
  }
  
  if(species=="mouse"){
    library(biomaRt)
    human=readRDS("/mnt/ssd1/longjie/newest_scripts/database/human_gene_ensembl.rds")
    mouse=readRDS("/mnt/ssd1/longjie/newest_scripts/database/mouse_gene_ensembl.rds")
    
    gene_mouse=row.names(total_matrix)
    print("Converting mouse gene symbols to human...")
    gene_human=getLDS(attributes = c("mgi_symbol"),
                      filters = "mgi_symbol",
                      values=gene_mouse,mart = mouse,
                      attributesL = c("hgnc_symbol"),
                      martL = human,uniqueRows = T)
    gene_human=gene_human[!duplicated(gene_human$MGI.symbol),]
    gene_human=gene_human[!duplicated(gene_human$HGNC.symbol),]
    total_matrix=total_matrix[gene_human$MGI.symbol,]
    row.names(total_matrix)=gene_human$HGNC.symbol
  }
  
  geneset=getGmt(gspath)
  print("step2. run gsva...")
  dim(total_matrix)
  total_es <- gsva(as.matrix(total_matrix), geneset, method=method,
                   min.sz=min.sz, max.sz=2000, verbose=FALSE, parallel.sz=parallel.sz)
  dim(total_es)
  ##find markers part----
  print("step3. find gs markers of clusters...")
  marker_gs=data.frame()
  for (i in levels(seurat_object)) {
    print(paste0("####-------",i," cluster------#####"))
    group1=i
    group_cell=supercell_cluster
    group_cell=ifelse(supercell_cluster==i,"group1","group2")
    group_cell=factor(group_cell,levels = c("group1","group2"))
    adjPvalueCutoff <- adjPvalueCutoff
    logFCcutoff <- logFCcutoff
    design <- model.matrix(~ factor(group_cell))
    vsname=paste0("group2","vs","group1")
    colnames(design) <- c("group1",vsname)
    fit <- lmFit(total_es, design)
    fit <- eBayes(fit)
    allGeneSets <- topTable(fit, coef=vsname,number=Inf)
    DEgeneSets <- topTable(fit, coef=vsname,number=Inf,p.value=adjPvalueCutoff, lfc = logFCcutoff,adjust="BH")
    #head(DEgeneSets)
    if(nrow(DEgeneSets)==0){
      stop("No geneset through the filter!")
    }
    DEgeneSets=subset(DEgeneSets,subset=logFC<0)
    DEgeneSets=DEgeneSets[order(DEgeneSets$t,decreasing = F),]
    count_up_gs=length(which(DEgeneSets$t<0))
    
    if(count_up_gs>=topn){
      DEgeneSets_top_name=head(row.names(DEgeneSets),topn)
    }else{
      DEgeneSets_top_name=head(row.names(DEgeneSets),count_up_gs)
    }
    marker_gs_i=data.frame(cluster=c(rep(i,length(DEgeneSets_top_name))),geneset=DEgeneSets_top_name)
    marker_gs=rbind(marker_gs,marker_gs_i)
  }
  marker_gs=subset(marker_gs,subset=!duplicated(marker_gs$geneset))
  write.csv(marker_gs,paste0("marker_geneset_",filename,".csv"))
  
  ##pheatmap part----
  escore=total_es
  if(!is.null(pathname)){
    escore=escore[pathname,]
  }else{
    escore=escore[marker_gs$geneset,]
  }
  escore_t=as.data.frame(t(escore))
  pheat_matrix=aggregate(escore_t,list(supercell_cluster),mean)
  Group.1=pheat_matrix$Group.1
  pheat_matrix=subset(pheat_matrix,select=-Group.1)
  escore_scale=as.data.frame(scale(pheat_matrix))
  escore_scale=as.data.frame(t(escore_scale))
  colnames(escore_scale)=Group.1
  escore_scale=escore_scale[,levels(seurat_object)]
  
  pdf(paste0("GSVA_cluster_diff_",filename,".pdf"),width = 1.15*width,height = height)
  p=pheatmap(escore_scale,fontsize_row = fontsize_row,cluster_rows = cluster_rows,cluster_cols = cluster_cols,
             treeheight_row = 0,
             color = col)
  print(p)
  dev.off()
  return(escore_scale)
}





