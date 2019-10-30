###Ananlyse scMatch output
library(ggplot2)

#Path to output from cellranger
cell_ranger_path

data <- Read10X(data.dir = paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix"))
data2 <- CreateSeuratObject(counts = data, project = "PROJECT_NAME", min.cells = 3, min.features = 200)



annot <- read.csv(paste0(cell_ranger_path, "/outs/filtered_feature_bc_matrix/annotation_result_keep_all_genes/mouse_Spearman_top_ann.csv"))

hist(annot$top.correlation.score, breaks=100, main="Vehicle")
cutoff <- summary(annot$top.correlation.score)[2]
abline(v=cutoff, col="red")

annot <- annot[annot$top.correlation.score > cutoff,]
annot$cell.type <- as.character(annot$cell.type)


df <- table(annot$cell.type)
df <- data.frame(Cell_type = names(df),
                       Number = as.vector(df))
df$Experiment <- "Vehicle"
total_cells <- sum(df$Number)
df$Percentage <- (df$Number/total_cells)*100


ggplot(df, aes(x=Cell_type, y=Percentage))+
  geom_bar(stat='identity', aes(fill=Experiment), position='dodge')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
