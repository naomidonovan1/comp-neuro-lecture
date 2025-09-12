###############################################
# DEMO: Intro to scRNA-seq 
# Part 1: R basics
# Part 2: Exploring scRNA-seq data
# Part 3: Analyzing scRNA-seq data
###############################################

# --------------------------
# ==== Part 1: R basics ====
# --------------------------

# Assign a number
x <- 5   # "<-" is preferred for assignment; "=" is mostly for function arguments
x

# Assign a vector
y <- c(1, 2, 3, 4, 5)
y

# Arithmetic
x * y
mean(y)

# Data types
num  <- 3.14         # numeric
txt  <- "astrocyte"  # character
bool <- TRUE         # logical
class(num); class(txt); class(bool)

# Looping
for (i in 1:5) {
  print(i)
}

# Data frames
df <- data.frame(
  cell   = c("cell1", "cell2", "cell3"),
  type   = c("Neuron", "Astrocyte", "Microglia"),
  counts = c(100, 200, 150)
)
df

# Subsetting
df$cell                       # one column
df[1, ]                       # first row
df[df$type == "Neuron", ]     # filter rows

# Built-in functions
sum(df$counts)
mean(df$counts)

# Define your own function
double <- function(x) { x * 2 }
double(10)

# Plotting with ggplot2
library(ggplot2)             
ggplot(df, aes(x = type, y = counts, fill = type)) +
  geom_col() +
  theme_classic() +
  ggtitle("Cell counts")
# ?ggplot  # check the help page for a function


# ----------------------------------------------------
# ==== Part 2: Exploring single-cell RNA-seq data ====
# ----------------------------------------------------

library(Matrix)
library(dplyr)
library(ggrepel)
library(Seurat)
library(SingleCellExperiment)
library(ExperimentHub)

# Create connection to Bioconductor ExperimentHub (public dataset repository)
eh <- ExperimentHub()
query(eh, "Zeisel")   # search for datasets with 'Zeisel' in metadata

# Load Zeisel 2015 mouse cortex+hippocampus dataset (counts + metadata)
counts   <- eh[["EH2580"]]
metadata <- eh[["EH2582"]]

# Terminology:
# rows = genes (e.g., Gfap, Gad1, Sst, Pvalb)
# columns = cells (samples)
# entries = counts (# of unique UMIs/reads from that gene in that cell)

# Convert to sparse matrix (saves memory, since scRNA-seq data are mostly zeros)
if (!inherits(counts, "dgCMatrix")) 
  counts <- Matrix(counts, sparse = TRUE)

# Align metadata rows to matrix columns
stopifnot(ncol(counts) == nrow(metadata))
if (!length(intersect(colnames(counts), rownames(metadata)))) 
  rownames(metadata) <- colnames(counts)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(
  assays  = list(counts = counts), 
  colData = metadata)
sce

# Convert to Seurat object
obj <- CreateSeuratObject(
  counts    = assay(sce, "counts"),
  meta.data = as.data.frame(colData(sce)))

# Standard preprocessing workflow
obj <- NormalizeData(obj) |>
  FindVariableFeatures(nfeatures = 2000) |>
  ScaleData() |>
  RunPCA(npcs = 30)

# Show variance explained (helps choose how many PCs to use)
ElbowPlot(obj, ndims = 30)

# Explore metadata categories
colnames(obj@meta.data)
table(obj$level1class)

# PCA scatterplot
DimPlot(obj, reduction = "pca", group.by = "level1class") +
  ggtitle("PCA: broad cell classes")

# UMAP + clustering (neighbors → UMAP → clusters)
obj <- obj |>
  FindNeighbors(dims = 1:30) |>
  RunUMAP(dims = 1:30) |>
  FindClusters(resolution = 0.4)

# Plot UMAPs
DimPlot(obj, reduction = "umap", group.by = "level1class", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: major cell types (level1class)")
DimPlot(obj, reduction = "umap", group.by = "level2class", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: cell subtypes (level2class)")
DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Seurat clusters")


# ----------------------------------------------------
# ==== Part 3: Analyzing single-cell RNA-seq data ====
# ----------------------------------------------------

# Marker gene expression
markers <- c("Slc17a7",  # excitatory neurons
             "Gad1",     # inhibitory neurons
             "Gfap",     # astrocytes
             "Mbp",      # oligodendrocytes
             "Cx3cr1")   # microglia

FeaturePlot(obj, features = markers[1]) + ggtitle(paste0(markers[1], " (Excitatory neurons)"))
FeaturePlot(obj, features = markers[2]) + ggtitle(paste0(markers[2], " (Inhibitory neurons)"))
FeaturePlot(obj, features = markers[3]) + ggtitle(paste0(markers[3], " (Astrocytes)"))
FeaturePlot(obj, features = markers[4]) + ggtitle(paste0(markers[4], " (Oligodendrocytes)"))
FeaturePlot(obj, features = markers[5]) + ggtitle(paste0(markers[5], " (Microglia)"))

# Cell type proportions
obj@meta.data |>
  dplyr::count(level1class, name = "n") |>
  dplyr::mutate(prop = n / sum(n)) |>
  ggplot(aes(x = reorder(level1class, prop), y = prop, fill = level1class)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Proportion of cells", 
       title = "Cell-type composition (level1class)") +
  theme_classic() + theme(legend.position = "none")

# Differential expression: pyramidal CA1 vs interneurons
Idents(obj) <- obj$level1class
de_test <- FindMarkers(obj, ident.1 = "pyramidal CA1", ident.2 = "interneurons",
                          test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1)
# de_test <- FindMarkers(obj, ident.1 = "Astro1", ident.2 = "Astro2",
#                           test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1)
de_test$gene <- rownames(de_test)

# Peek at top results
head(de_test[order(de_test$p_val_adj), ], 10)

# Volcano plot
lfc_cut  <- 0.5
padj_cut <- 0.05
de_test$neglog10_padj <- -log10(de_test$p_val_adj + 1e-300)
de_test$sig <- de_test$p_val_adj < padj_cut & abs(de_test$avg_log2FC) >= lfc_cut

# Select top genes for labeling
top_pos <- de_test %>% 
  filter(avg_log2FC > 0) %>% 
  arrange(desc(neglog10_padj)) %>% 
  slice_head(n = 5)
top_neg <- de_test %>% 
  filter(avg_log2FC < 0) %>% 
  arrange(desc(neglog10_padj)) %>% 
  slice_head(n = 5)
top_genes <- bind_rows(top_pos, top_neg)

# Plot volcano
ggplot(de_test, aes(x = avg_log2FC, y = neglog10_padj)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dotted") +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = Inf) +
  labs(x = "log2 FC",
       y = "-log10(adj p)", title = "DEGs: pyramidal CA1 vs interneurons") +
  theme_classic()