###
library(clusterProfiler)
library(org.Hs.eg.db)

CD4_degs <- FindMarkers(CD4.combined,ident.1 = 'CD4_Tcm_TOX',
            ident.2 = 'CD4_Tcm_PASK',group.by = 'Celltype3',
            assay = 'RNA',logfc.threshold = 0.25,min.pct = 0.1)
CD4_degs$Gene <- rownames(CD4_degs)


ids=bitr(rownames(subset(CD4_degs,avg_log2FC>0)),'SYMBOL','ENTREZID','org.Hs.eg.db') 
ekegg <- enrichKEGG(gene = ids$ENTREZID,
                    organism = "hsa")
#dotplot(ekegg,showCategory = 10,font.size = 10,title = "KEGG up with age")
ekegg_df_up <- data.frame(ekegg)
ekegg_df_up$Group <- 'CD4_Tcm_TOX'
#
ids=bitr(rownames(subset(CD4_degs,avg_log2FC<0)),'SYMBOL','ENTREZID','org.Hs.eg.db') 
ekegg <- enrichKEGG(gene = ids$ENTREZID,
                    organism = "hsa")
#dotplot(ekegg,showCategory = 10,font.size = 10,title = "KEGG down with ag")
ekegg_df_down <- data.frame(ekegg)
ekegg_df_down$Group <- 'CD4_Tcm_PASK'

data <- rbind(ekegg_df_up,ekegg_df_down) %>% 
  mutate(p2 = ifelse(Group == "CD4_Tcm_TOX", -log10(p.adjust), log10(p.adjust))) %>% 
  arrange(Group,p2) 

paths <- c(head(data,6)$Description,tail(data,6)$Description)

###
library(GSVA)
library(limma)

de_gsva  <- function(exprSet,meta,compare = NULL){
  
  
  allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  
  fit <- lmFit(exprSet,design)
  if(length(unique(meta))==2){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  
  return(allDiff)
}

Group_GSVA <- function(seurat_obj, Group, category='KEGG',compare = compare) {
  expr=as.matrix(seurat_obj@assays$RNA@data)
  
  if(category == 'GO.BP'){                          
    msgdC5 = msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
    GO.BPSet = msgdC5 %>% split(x = .$gene_symbol, f = .$gs_name)
    sc_GSVA <- gsva(expr, gset.idx.list = GO.BPSet, kcdf="Gaussian",method = "gsva",
                    parallel.sz=30)   
  }
  if(category == 'KEGG'){
    kegg_df <- read.csv("KEGGREST_WithGene.csv",row.names = 1)
    gene_list <- strsplit(kegg_df$hgnc_symbol,split = ',')
    names(gene_list) <- kegg_df$pathway_name
    sc_GSVA <- gsva(expr, gset.idx.list =gene_list, kcdf="Gaussian",method = "gsva",
                    parallel.sz=30)   
  }
  
  
  #
  meta <- seurat_obj@meta.data[,Group]
  meta <- str_replace(meta,pattern = '-',replacement = '_')
  Diff =de_gsva(exprSet = sc_GSVA,meta = meta,compare = compare)
  Padj_threshold=0.05
  idiff <-Diff[[compare]]
  df <- data.frame(ID = rownames(idiff), score = idiff$t )
  df$group =sapply(1:nrow(idiff),function(x){
    if(idiff[x,"logFC"]>0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("up")}
    else if(idiff[x,"logFC"]<0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("down")}
    else{return("noSig")}
  })
  
  df$hjust = ifelse(df$score>0,1,0)
  df$nudge_y = ifelse(df$score>0,-0.1,0.1)
  sortdf <- df[order(df$score),]
  sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
  limt = max(abs(df$score))
  Mono_sortdf_sub <- sortdf
  p <- ggplot(Mono_sortdf_sub, aes(ID, score,fill=group)) + 
    geom_bar(stat = 'identity',alpha = 0.7) + 
    scale_fill_manual(breaks=c("down","noSig","up"),
                      values = c("#008020","grey","#08519C"))+
    geom_text(data = Mono_sortdf_sub, aes(label = Mono_sortdf_sub$ID, y = Mono_sortdf_sub$nudge_y),
              nudge_x =0,nudge_y =0,hjust =Mono_sortdf_sub$hjust,
              size = 3)+
    labs(x = paste0('KEGG'," pathways"),
         title = 'Monocyte')+
    scale_y_continuous(limits=c(-25,25))+
    coord_flip() + 
    theme_bw() + 
    theme(panel.grid =element_blank())+
    theme(panel.border = element_rect(size = 0.6)
          #panel.border = element_blank()
    )+
    theme(plot.title = element_text(hjust = 0.5,size = 18),
          axis.text.y = element_blank(),
          axis.title = element_text(hjust = 0.5,size = 18),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = limt
    )
  return(list(sortdf,p))
}
Mono_kegg1 <-Group_GSVA(mono_sub,Group = 'Group',category = "KEGG",compare = 'HC-HLH')
