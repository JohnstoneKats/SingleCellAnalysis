library(dplyr)
library(Seurat)
library(reticulate)
library(sctransform)
library(SeuratData)
library(cowplot)
library(ggplot2)
use_virtualenv(".virtualenvs/UMAP")

#Path to output from cellranger
cell_ranger_path_drug
cell_ranger_path_veh


drug.data <- Read10X(data.dir = paste0(cell_ranger_path_drug, "/outs/filtered_feature_bc_matrix"))
veh.data <- Read10X(data.dir = paste0(cell_ranger_path_ve, "/outs/filtered_feature_bc_matrix"))

drug <- CreateSeuratObject(counts = drug.data, project = "PROJECT_NAME", min.cells = 3, min.features = 200)
drug  #14419 features, 3950 samples
drug[["percent.mt"]] <- PercentageFeatureSet(drug, pattern = "^mt-")
VlnPlot(drug, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
drug <- subset(drug, subset = nFeature_RNA > 200 & percent.mt < 5)

veh <- CreateSeuratObject(counts = veh.data, project = "PROJECT_NAME", min.cells = 3, min.features = 200)
veh  #14028 features, 4077 samples
veh[["percent.mt"]] <- PercentageFeatureSet(veh, pattern = "^mt-")
VlnPlot(veh, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
veh <- subset(veh, subset = nFeature_RNA > 200 & percent.mt < 5)

# Identify the 10 most highly variable genes
top10_drug <- head(VariableFeatures(drug), 10)
top10_veh <- head(VariableFeatures(veh), 10)

# run sctransform
drug %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE) -> drug

DimPlot(label = TRUE) + NoLegend()


veh %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE) -> veh
  

#Load scMatch annotations ('annot' object) generated from the scMatch_analysis.R script

load("annot_veh.Rdata")
load("annot_drug.Rdata")

annot_veh$top.sample <- NULL
annot_veh$top.correlation.score <- NULL
annot_veh$cell <- gsub("-1", "", LL_38_annot$cell)

temp_names <- names(veh$SCT_snn_res.0.8)
veh$cell_types_scMatch <- LL_38_annot$cell.type[match(temp_names, annot_veh$cell)]

veh$cell_types_scMatch2 <- veh$cell_types_scMatch
veh$cell_types_scMatch2 <- gsub(", fetal liver derived", "", veh$cell_types_scMatch2)
veh$cell_types_scMatch2 <- gsub(", placenta derived", "", veh$cell_types_scMatch2)
veh$cell_types_scMatch2 <- gsub("- alternatively activated", "", veh$cell_types_scMatch2)

DimPlot(veh, label = TRUE, group.by = 'cell_types_scMatch2')
DimPlot(veh, label = FALSE, group.by = 'cell_types_scMatch2')


annot_drug$top.sample <- NULL
annot_drug$top.correlation.score <- NULL
annot_drug$cell <- gsub("-1", "", LL_30_annot$cell)

temp_names <- names(drug$SCT_snn_res.0.8)
drug$cell_types_scMatch <- LL_30_annot$cell.type[match(temp_names, annot_drug$cell)]

drug$cell_types_scMatch2 <- drug$cell_types_scMatch
drug$cell_types_scMatch2 <- gsub(", fetal liver derived", "", drug$cell_types_scMatch2)
drug$cell_types_scMatch2 <- gsub(", placenta derived", "", drug$cell_types_scMatch2)
drug$cell_types_scMatch2 <- gsub("- alternatively activated", "", drug$cell_types_scMatch2)

DimPlot(drug, label = TRUE, group.by = 'cell_types_scMatch2')
DimPlot(drug, label = FALSE, group.by = 'cell_types_scMatch2')


###Combining plots
combined.integrated <- RunPCA(combined.integrated, verbose = FALSE)
combined.integrated <- RunUMAP(combined.integrated, dims = 1:30)
DimPlot(veh, label = TRUE)

options(future.globals.maxSize = 4000 * 1024^2)
drug <- CreateSeuratObject(counts = drug.data, project = "drug", min.cells = 3, min.features = 200)
drug
drug[["percent.mt"]] <- PercentageFeatureSet(drug, pattern = "^mt-")
VlnPlot(drug, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
drug <- subset(drug, subset = nFeature_RNA > 200 & percent.mt < 5)

veh <- CreateSeuratObject(counts = veh.data, project = "veh", min.cells = 3, min.features = 200)
veh
veh[["percent.mt"]] <- PercentageFeatureSet(veh, pattern = "^mt-")
VlnPlot(veh, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
veh <- subset(veh, subset = nFeature_RNA > 200 & percent.mt < 5)



##Combination
combined <- merge(veh, y = drug, add.cell.ids = c("veh", "drug"), project = "PROJECT_NAME")

combined.list <- SplitObject(combined)

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE)
}

genes_veh <- rownames(veh)
genes_drug <- rownames(drug)
genes_veh[genes_veh %in% genes_drug]


combined.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 3000)
combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = combined.features, 
                                    verbose = FALSE)
combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                           anchor.features = combined.features, verbose = FALSE)
combined.integrated <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
combined.integrated <- RunPCA(combined.integrated, verbose = FALSE)
combined.integrated <- RunUMAP(combined.integrated, dims = 1:30)
plots <- DimPlot(combined.integrated, combine = FALSE)


plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
                                                                                                              byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)

veh_cells <- WhichCells(object=combined.integrated, ident="veh")
drug_cells <- WhichCells(object=combined.integrated, ident="drug")

combined.integrated %>% 
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE, resolution=1.1) -> combined.integrated2 #1.1
  #DimPlot(label=TRUE)

DimPlot(combined.integrated2, label = TRUE)

markers <- FindAllMarkers(combined.integrated2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


##Adding the annotations from scMatch to the combined plot
DimPlot(combined.integrated2, reduction = "umap", split.by = "orig.ident")

scMatch_annot <- rbind(annot_drug, annot_veh)

combined.integrated2$integrated_snn_res.1.1
temp_names <- names(combined.integrated2$cell_types_scMatch)
combined.integrated2$cell_types_scMatch <- scMatch_annot$cell.type[match(temp_names, scMatch_annot$cell)]

combined.integrated2$cell_types_scMatch2 <- combined.integrated2$cell_types_scMatch
combined.integrated2$cell_types_scMatch2 <- gsub(", fetal liver derived", "", combined.integrated2$cell_types_scMatch2)
combined.integrated2$cell_types_scMatch2 <- gsub(", placenta derived", "", combined.integrated2$cell_types_scMatch2)
combined.integrated2$cell_types_scMatch2 <- gsub(" - alternatively activated", "", combined.integrated2$cell_types_scMatch2)

combined.integrated2$orig.ident <- factor(combined.integrated2$orig.ident, levels=c("veh", "drug"))

DimPlot(combined.integrated2, reduction = "umap", group.by = "cell_types_scMatch2",
        split.by='orig.ident')

DimPlot(combined.integrated2, reduction = "umap", group.by = "orig.ident",
        split.by = "cell_types_scMatch2")

which_no_na <- names(combined.integrated2$cell_types_scMatch2)[!is.na(combined.integrated2$cell_types_scMatch2)]
combined.integrated_no_na <- subset(combined.integrated2, cells=which_no_na)

DimPlot(combined.integrated_no_na, reduction = "umap", group.by = "orig.ident",
        split.by = "cell_types_scMatch2")

#Calculation of cell proportions
cell_proportions <- vector()
for(cluster in 0:14){
  cluster_cells <- WhichCells(object=combined.integrated2, ident=cluster)
  n_veh <- sum(cluster_cells %in% veh_cells)
  n_drug <- sum(cluster_cells %in% drug_cells)
  
  p_veh <- (n_veh/length(cluster_cells))*100
  p_drug <- (n_drug/length(cluster_cells))*100
  
  cell_proportions <- rbind(cell_proportions,
                            c(cluster, length(cluster_cells), n_veh, n_drug, p_veh, p_drug))
}

colnames(cell_proportions) <- c("Cluster",
                                "Total_cells", "N_veh", "N_drug", "P_veh", "P_drug")



cell_proportions <- as.data.frame(cell_proportions)
significant_genes$P_veh <- cell_proportions$P_veh[match(significant_genes$cluster, cell_proportions$Cluster)]
significant_genes$P_drug <- cell_proportions$P_drug[match(significant_genes$cluster, cell_proportions$Cluster)]
