library(destiny)
library(Seurat)
library(Biobase)
library(readxl)

#https://bioconductor.org/packages/release/bioc/html/destiny.html

#Path to output from cellranger
cell_ranger_path

data <- Read10X(data.dir = paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix"))
data2 <- CreateSeuratObject(counts = data, project = "PROJECT_NAME", min.cells = 3, min.features = 200)

barcodes <- read.delim(paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix/barcodes.tsv"), header=FALSE)
features <- read.delim(paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix/features.tsv"), header=FALSE)

matrix_data <- as.matrix(GetAssayData(data2, slot = "counts"))

ct <- as.ExpressionSet(as.data.frame(matrix_data))

dm <- DiffusionMap(ct)
plot(dm)

