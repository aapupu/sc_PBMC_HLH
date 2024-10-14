library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
data.input =GetAssayData(object =pbmc.Non_doublet,slot = 'data')# normalized data matrix
meta = pbmc.Non_doublet@meta.data # a dataframe with rownames containing cell mata data

cell.use_HC = rownames(meta)[meta$Group  == "HC"] # extract the cell names from disease data
cell.use_HLH = rownames(meta)[meta$Group  == "HLH"] # extract the cell names from disease data

##
meta_HC = meta[cell.use_HC, ]
meta_HC$Celltype3 = droplevels(meta_HC$Celltype3, exclude = setdiff(levels(meta_HC$Celltype3),unique(meta_HC$Celltype3)))

cellchat_HC <- createCellChat(object = data.input[, cell.use_HC], meta = meta_HC, group.by = "Celltype3")
cellchat_HLH <- createCellChat(object =  data.input[, cell.use_HLH], meta = meta[cell.use_HLH, ], group.by = "Celltype3")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling",'Cell-Cell Contact','ECM-Receptor'))
cellchat_HC@DB <- CellChatDB.use
cellchat_HLH@DB <- CellChatDB.use

cellchat_HC <- subsetData(cellchat_HC) 
cellchat_HLH <- subsetData(cellchat_HLH) 

future::plan("multiprocess", workers = 30)
#
cellchat_HC <- identifyOverExpressedGenes(cellchat_HC)
cellchat_HC <- identifyOverExpressedInteractions(cellchat_HC)

cellchat_HLH <- identifyOverExpressedGenes(cellchat_HLH)
cellchat_HLH <- identifyOverExpressedInteractions(cellchat_HLH)
#
cellchat_HC <- computeCommunProb(cellchat_HC)
cellchat_HC <- filterCommunication(cellchat_HC, min.cells = 10)

cellchat_HLH <- computeCommunProb(cellchat_HLH)
cellchat_HLH <- filterCommunication(cellchat_HLH, min.cells = 10)
#
cellchat_HC <- computeCommunProbPathway(cellchat_HC)
cellchat_HLH <- computeCommunProbPathway(cellchat_HLH)
#
cellchat_HC <- aggregateNet(cellchat_HC)
cellchat_HLH <- aggregateNet(cellchat_HLH)

cellchat_HC <- netAnalysis_computeCentrality(cellchat_HC, slot.name = "netP")
cellchat_HLH <- netAnalysis_computeCentrality(cellchat_HLH, slot.name = "netP")
