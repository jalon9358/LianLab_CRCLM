#### mycellchat 20211102 by longjie ####
library(future.apply)
library(CellChat)
library(patchwork)

###function part----
mycellchat=function(object,database=CellChatDB.human,db.use="Secreted Signaling"){
  data.input=object@assays$RNA@data
  meta=data.frame(labels=object@active.ident,row.names = colnames(object))

  cellchat=createCellChat(object = data.input,meta = meta,group.by = "labels")
  cellchatdb=database
  # use a subset of CellChatDB for cell-cell communication analysis
  cellchatdb.use <- subsetDB(cellchatdb, search = db.use)

  #cellchatdb.use$interaction=cellchatdb.use$interaction[-c(717,719,720),]
  
  cellchat@DB=cellchatdb.use
  cellchat=subsetData(cellchat)
  plan(strategy = "multisession", workers = 6) ## Run in parallel on local computer
  options(future.globals.maxSize = 36*1024^3)
  cat("It will spend some times to run, please wait...\n")
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  plan(strategy = "multisession", workers = 1)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  result=list(df.net=df.net,cellchat=cellchat)
  return(result)
}



