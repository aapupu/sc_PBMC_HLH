library(nichenetr)


ligand_target_matrix = readRDS("/nichenet/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]# target genes in rows, ligands in columns

lr_network = readRDS("/nichenet/lr_network.rds")
head(lr_network)

weighted_networks = readRDS("/nichenet/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig)# interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
#########
DE_table_receiver <- subset(mono_Group,avg_log2FC>0) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

Myeloid.combined2@active.ident <- as.factor(Myeloid.combined2$Group)
expressed_genes_receiver = get_expressed_genes('HLH', Myeloid.combined2, pct = 0.10,assay_oi = 'RNA')
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# geneset_oi = rownames(subset(cDC_degs[cDC_degs$cluster %in% 'migDC_LAMP3',],p_val_adj<0.05 & avg_log2FC>0.25)) %>% .[. %in% rownames(ligand_target_matrix)]


expressed_ligands = ligands
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(dplyr::desc(aupr_corrected)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
# best_upstream_ligands <- c(best_upstream_ligands,'IFNG')
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes of Monocyte (HLH)", color = "#B2182B",
                                          legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic"),axis.text=element_text(color='black')) +
  scale_fill_gradient2(low = "whitesmoke",  high = "#B2182B", breaks = c(0,0.0070))
p_ligand_target_network
