# ==========================================
# 1. DATA EXPLORATION AND QUALITY CONTROL
# ==========================================

# Load packages
library(ggplot2)
library(edgeR)

# Load metadata and gene expression counts
metadata <- read.csv("data/metadata.csv")

counts <- read.csv(
  "data/counts.csv",
  row.names = 1,
  check.names = FALSE
)

# Check dimensions
dim(metadata)
dim(counts)

# Check data types
str(metadata)

# Check whether sample names match
all(metadata$description == colnames(counts))

# Check missing values
colSums(is.na(metadata))

# Check unique donors
length(unique(metadata$donor_id))

# Check diagnosis distribution
table(metadata$diagnosis)

# Check brain region distribution
table(metadata$brain_region)

# Check sex distribution
table(metadata$sex)

# Check age, PMI and brain pH
summary(metadata$age)
summary(metadata$PMI)
summary(metadata$brain_pH)# ==========================================
# 2. LIBRARY SIZE AND GENE EXPRESSION QC
# ==========================================

# Calculate library sizes
library_size <- colSums(counts)

# Summary of library sizes
summary(library_size)

# Plot library sizes
boxplot(
  library_size / 1e6,
  main = "Library sizes per sample",
  ylab = "Millions of reads",
  col = "lightblue"
)

# Distribution of total gene counts
hist(
  log10(rowSums(counts) + 1),
  breaks = 50,
  main = "Distribution of total gene counts",
  xlab = "Log10(total counts + 1)",
  col = "steelblue"
)

# Count genes with zero expression
sum(rowSums(counts) == 0)# ==========================================
# 3. FILTERING LOW-EXPRESSED GENES
# ==========================================

# Create edgeR object
dge <- edgeR::DGEList(counts = as.matrix(counts))

# Define groups based on brain region and diagnosis
group <- interaction(
  metadata$brain_region,
  metadata$diagnosis
)

# Filter low-expressed genes using CPM
keep <- edgeR::filterByExpr(dge, group = group)

# Keep filtered genes
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

# Check number of genes after filtering
dim(dge_filtered)

# Calculate number of removed genes
sum(!keep)# ==========================================
# 4. TMM NORMALIZATION AND LOG TRANSFORMATION
# ==========================================

# Perform TMM normalization
dge_filtered <- edgeR::normLibSizes(
  dge_filtered,
  method = "TMM"
)

# Inspect normalization factors
summary(dge_filtered$samples$norm.factors)

# Transform counts to log2-CPM
logCPM <- edgeR::cpm(
  dge_filtered,
  log = TRUE,
  prior.count = 2
)

# Check dimensions
dim(logCPM)# ==========================================
# 5. PCA ANALYSIS
# ==========================================

# Perform PCA
pca <- prcomp(t(logCPM), scale. = FALSE)

# Calculate explained variance
variance_explained <- pca$sdev^2 / sum(pca$sdev^2)

# Display variance explained by first five PCs
round(variance_explained[1:5] * 100, 2)

# Prepare PCA data
pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  diagnosis = metadata$diagnosis,
  brain_region = metadata$brain_region,
  age = metadata$age,
  PMI = metadata$PMI
)

# PCA colored by diagnosis
ggplot(pca_data, aes(x = PC1, y = PC2, color = diagnosis)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA of gene expression by diagnosis",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

# PCA colored by brain region
ggplot(pca_data, aes(x = PC1, y = PC2, color = brain_region)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA of gene expression by brain region",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

# PCA colored by age
ggplot(pca_data, aes(x = PC1, y = PC2, color = age)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal()

# PCA colored by PMI
ggplot(pca_data, aes(x = PC1, y = PC2, color = PMI)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal()# ==========================================
# 6. SAVE RESULTS
# ==========================================

# Create results folder if it does not exist
dir.create("results", showWarnings = FALSE)

# Save normalized expression data
write.csv(
  logCPM,
  "results/logCPM_normalized.csv"
)

# Save filtered and normalized edgeR object
saveRDS(
  dge_filtered,
  "results/edgeR_filtered_normalized.rds"
)

# Save PCA results
saveRDS(
  pca,
  "results/PCA.rds"
)