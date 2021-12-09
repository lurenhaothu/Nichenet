library(nichenetr)
library(tidyverse)

#Step 0: Load required packages, NicheNet’s ligand-target prior model 
#and processed expression data of interacting cells
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#ligand_target_matrix = read.csv("ligand_target_matrix.csv")
#ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

TM_expression = read.csv("TM expression.csv")
SC_expression = read.csv("SC expression.csv")

### Step 1: Define expressed genes in sender and receiver cell populations

expressed_genes_sender = TM_expression$ï..Gene
expressed_genes_receiver = SC_expression$Gene

#length(expressed_genes_sender)
#length(expressed_genes_receiver)

###Step 2: Define the gene set of interest and a background of genes

# reference:: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3475402/

geneset_oi = read.csv("gene_oi_large.csv")$gene %>% .[. %in% rownames(ligand_target_matrix)] %>% .[. %in% SC_expression$Gene]
# head(geneset_oi)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
#head(background_expressed_genes)

###Step 3: Define a set of potential ligands

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

ligands = lr_network %>% pull(from) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()

expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) 
# head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
# head(potential_ligands)

###Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities %>% arrange(-pearson) 

best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(15, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

###Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df[-27,]

#nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

# nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot(
  "Prioritized TM-ligands","Junction genes in SC cells", 
  color = "purple",legend_position = "top", 
  x_axis_position = "top",legend_title = "Regulatory potential") + 
  #scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + 
  theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network

###Follow-up analysis 1: Ligand-receptor network inference for top-ranked ligands

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized TM-ligands","Receptors expressed by SC endothelial cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

###Follow-up analysis 2: Visualize expression of top-predicted ligands and their target genes in a combined heatmap

library(RColorBrewer)
library(cowplot)
library(ggpubr)

ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
p_ligand_pearson

# expression_df_CAF = expression[CAF_ids,order_ligands] %>% 
#   data.frame() %>% 
#   rownames_to_column("cell") %>% 
#   as_tibble() %>% 
#   inner_join(sample_info %>% select(cell,tumor), by =  "cell")
# 
# aggregated_expression_CAF = expression_df_CAF %>% group_by(tumor) %>% select(-cell) %>% summarise_all(mean)
# 
# aggregated_expression_df_CAF = aggregated_expression_CAF %>% select(-tumor) %>% t() %>% magrittr::set_colnames(aggregated_expression_CAF$tumor) %>% data.frame() %>% rownames_to_column("ligand") %>% as_tibble() 
# 
# aggregated_expression_matrix_CAF = aggregated_expression_df_CAF %>% select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_CAF$ligand)
# 
# order_tumors = c("HN6","HN20","HN26","HN28","HN22","HN25","HN5","HN18","HN17","HN16") # this order was determined based on the paper from Puram et al. Tumors are ordered according to p-EMT score.

aggregated_expression_matrix_TM = TM_expression %>% select(Value) %>% as.matrix() %>% magrittr::set_rownames(TM_expression$ï..Gene)
aggregated_expression_matrix_TM_1 = aggregated_expression_matrix_TM[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Expression")

vis_ligand_tumor_expression = aggregated_expression_matrix_TM_1

#p_ligand_pearson = vis_ligand_tumor_expression %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
#p_ligand_pearson

library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
p_ligand_tumor_expression = vis_ligand_tumor_expression %>% make_heatmap_ggplot("Prioritized TM-ligands","Ligand expression", color = color[100],legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(averaged over\nsingle cells)") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_tumor_expression

# expression_df_target = expression[malignant_ids,geneset_oi] %>% data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% inner_join(sample_info %>% select(cell,tumor), by =  "cell") 
# 
# aggregated_expression_target = expression_df_target %>% group_by(tumor) %>% select(-cell) %>% summarise_all(mean)
# 
# aggregated_expression_df_target = aggregated_expression_target %>% select(-tumor) %>% t() %>% magrittr::set_colnames(aggregated_expression_target$tumor) %>% data.frame() %>% rownames_to_column("target") %>% as_tibble() 
# 
# aggregated_expression_matrix_target = aggregated_expression_df_target %>% select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)
# 
# vis_target_tumor_expression_scaled = aggregated_expression_matrix_target %>% t() %>% scale_quantile() %>% .[order_tumors,order_targets]

aggregated_expression_matrix_target = SC_expression %>% select(Value) %>% as.matrix() %>% magrittr::set_rownames(SC_expression$Gene)
aggregated_expression_matrix_target_1 = aggregated_expression_matrix_target[order_targets,] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Expression")

vis_target_tumor_expression_scaled = t(aggregated_expression_matrix_target_1)

p_target_tumor_scaled_expression = vis_target_tumor_expression_scaled  %>% make_threecolor_heatmap_ggplot("Target expression","Target", low_color = color[1],mid_color = color[50], mid = 0.5, high_color = color[100], legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))
p_target_tumor_scaled_expression


figures_without_legend = plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  p_ligand_tumor_expression + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""), 
  NULL,
  NULL,
  p_target_tumor_scaled_expression + theme(legend.position = "none", axis.ticks = element_blank()) + xlab(""), 
  align = "hv",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_pearson)+ 6, ncol(vis_ligand_tumor_expression)+6, ncol(vis_ligand_target)) -2,
  rel_heights = c(nrow(vis_ligand_pearson), nrow(vis_target_tumor_expression_scaled) + 3)) 

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_pearson)),
  as_ggplot(get_legend(p_ligand_tumor_expression)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_target_tumor_scaled_expression)),
  nrow = 2,
  align = "h")

plot_grid(figures_without_legend, 
          legends, 
          rel_heights = c(10,2), nrow = 2, align = "hv")

vis_target_tumor_expression_scaled

