library(reticulate)
library(monocle3)
use_virtualenv(".virtualenvs/UMAP")

#Path to output from cellranger
cell_ranger_path

cds <- load_cellranger_data(cell_ranger_path)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

#UMAP
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=c(10^seq(-6,-1)))

plot_cells(cds)


#tSNE
cds_tSNE <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds_tSNE, reduction_method="tSNE")

cds <- learn_graph(cds)
plot_cells(cds,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
