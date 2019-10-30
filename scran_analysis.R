##scran test

library(scran)
library(dplyr)
library(Seurat)
library(reticulate)
library(sctransform)
library("org.Mm.eg.db")

cell_ranger_path

data <- Read10X(data.dir = paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix"))
data2 <- CreateSeuratObject(counts = data, project = "PROJECT_NAME", min.cells = 3, min.features = 200)

barcodes <- read.delim(paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix/barcodes.tsv"), header=FALSE)
features <- read.delim(paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix/features.tsv"), header=FALSE)

matrix_data <- as.matrix(GetAssayData(data2, slot = "counts"))


ensembl <- mapIds(org.Mm.eg.db, keys = rownames(matrix_data), column = "ENSEMBL", keytype="SYMBOL")
rownames(matrix_data) <- unname(ensembl)

sce <- SingleCellExperiment(list(counts=matrix_data))

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assigned <- cyclone(sce, pairs=mm.pairs)

head(assigned$scores)

table(assigned$phases)

