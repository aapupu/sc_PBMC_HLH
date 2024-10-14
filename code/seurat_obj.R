library(Seurat)
pbmc = CreateSeuratObject(counts,  project = "HLH", min.cells = 10, min.features = 200)
pbmc <- subset(pbmc, subset = percent.mt < 10)
pbmc_list <- SplitObject(pbmc,split.by = "orig.ident")
pbmc_list <- lapply(X = pbmc_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = pbmc_list)

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc_list, anchor.features = features)
# this command creates an 'integrated' data assay
pbmc.combined <- IntegrateData(anchorset = pbmc.anchors)

DefaultAssay(pbmc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE, vars.to.regress =c("nFeature_RNA", "nCount_RNA", "percent.mt"))

pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = FALSE)
DimPlot(pbmc.combined,group.by = 'orig.ident',reduction = 'pca')
ElbowPlot(pbmc.combined,ndims = 30)
#
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- RunTSNE(pbmc.combined, reduction = "pca", dims = 1:10)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 1.35) 
degs <- FindAllMarkers(pbmc.combined,assay = 'RNA',logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)

#
pbmc.combined$EBV <- ifelse(test = (colnames(pbmc.combined) %in% colnames(EBV_counts_sub)[colSums(EBV_counts_sub)!=0]),yes = 'Y',no = 'N')
DimPlot(pbmc.combined,cells.highlight = colnames(pbmc.combined)[pbmc.combined$EBV =='Y'],
        raster = F,sizes.highlight = 0.1,reduction = 'tsne')+ggtitle('EBV')+
  NoLegend()+theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
                   axis.text = element_blank(),
                   axis.ticks = element_blank())
