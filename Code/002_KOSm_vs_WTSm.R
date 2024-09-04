# The purpose of this script is to analyze the KO Large vs WT Large comparison.
# Note, that the actual running of DESeq2 was in the old Rproj "Delta_14".

# Clear environment
rm(list = ls())

# Load libraries
library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(magick)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Create output directories
opDir <- "Outputs/002_KOSm_vs_WTSm_Outputs"
if (!dir.exists(opDir)) {
  dir.create(opDir)
}

# Create subdirectory to hold rds files
rdsDir <- paste(opDir, "Rds_Files", sep = "/")
if (!dir.exists(rdsDir)) {
  dir.create(rdsDir)
}

# Create a subdirectory to hold DESeq2 results
ddsDir <- paste(opDir, "DESeq2", sep = "/")
if (!dir.exists(ddsDir)) {
  dir.create(ddsDir)
}

# Create a subdirectory for Figures
figDir <- paste(opDir, "Figures", sep = "/")
if (!dir.exists(figDir)) {
  dir.create(figDir)
}

#################################################################################
#----- Running Differential Expression Analysis

# Read in the counts
counts <- read.csv("Data_Files/Master_D14_Counts.csv")
rownames(counts) <- counts$X
counts$X <- NULL
metadata <- read.csv("Data_Files/Master_D14_Metadata.csv", header = TRUE, sep = ",")
rownames(metadata) <- metadata$X

# For this comparison, we are intersted in Groups C and D (KO Large and KO Small respectively)
meta <- metadata[metadata$Group == "D" | metadata$Group == "B",]

# Specify comarison name, file basenames, reference and treatment groups
comparison <- "KO Small vs WT Small"
fileBase <- "KO_Small_vs_WT_Small"
reference <- "WT_Small"
treatment <- "KO_Small"

# Filter the counts based on the keep vector
counts <- counts[,colnames(counts) %in% rownames(meta)]

# Check that everything is in the same order between counts and metadata
colnames(counts)[!(colnames(counts) %in% rownames(meta))]
all(colnames(counts) == rownames(meta))

# Prepare the DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Class)

# Set a reference level (KO Small)
dds$Class <- relevel(dds$Class, ref = reference)

# Pre-filter for low counts
smallestGroup <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroup
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Save rds file
saveRDS(dds, file = "Outputs/002_KOSm_vs_WTSm_Outputs/Rds_Files/KOSm_vs_WTSm_dds.Rds")


# Unhash the following line if you are reading in the dds file to edit figures
dds <- readRDS("Outputs/002_KOSm_vs_WTSm_Outputs/Rds_Files/KOSm_vs_WTSm_dds.Rds")

#################################################################################
#----- Principal Component Analysis

# Sanity check that the comparison was run correctly. Output should be "Class_WT_Large_vs_Class_WT_Small"
resultsNames(dds)

# Variance stabilize transform the data
vsd <- vst(dds)

# Run PCA and return data for custom plotting
PCA <- plotPCA(vsd, intgroup = "Class", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$group <- ifelse(PCA$group == "WT_Small", "WT", "KO")
PCA$group <- factor(PCA$group, levels = c("WT", "KO"))

# Save PCA as a csv
write.csv(PCA, file = "Outputs/002_KOSm_vs_WTSm_Outputs/DESeq2/KOSm_vs_WTSm_PrincipalComponents.csv")

# Set custom colors
colors <- c("WT" = "steelblue2", "KO" = "firebrick2")

# Plot PCA
plot <- ggplot(PCA, aes(PC1, PC2, fill = group, color = group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text(size = 4, aes(label = name), hjust = 1, vjust = 1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  coord_fixed() +
  theme_bw() +
  labs(title = "Small",
       color = "",
       fill = "") +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26, face = "bold"),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26, face = "bold"),
        legend.text = element_text(size = 24),
        title = element_text(size = 18))
ggsave("Outputs/002_KOSm_vs_WTSm_Outputs/Figures/KOSm_vs_WTSm_PCA_plot.tiff", plot, width = 10, height = 10, dpi = 300)

#################################################################################
#----- Extract differential expression analysis results
results <- as.data.frame(results(dds))
counts <- as.data.frame(counts(dds, normalized = TRUE))
results <- merge(results, counts, by = 0)
rownames(results) <- results$Row.names
results$Row.names <- NULL

# Get symbols
results$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(results))

# Remove unknown genes 
results <- results[!grepl("^ Gm\\d+$", results$Symbol),]
results <- results[!grepl("Rik$", results$Symbol),]
results <- results[!grepl("Rik\\d+$", results$Symbol),]
results <- results[!grepl("^ LOC", results$Symbol),]
results <- results[!grepl("AK\\d+$", results$Symbol),]
results <- results[!grepl("AY\\d+$", results$Symbol),]

# Remove IG genes
results <- results[!grepl("^ Igkv", results$Symbols),]
results <- results[!grepl("^ Ighg*", results$Symbols),]
results <- results[!grepl("^ Igha*", results$Symbols),]
results <- results[!grepl("^ Ighd*", results$Symbols),]
results <- results[!grepl("^ Ighe*", results$Symbols),]
results <- results[!grepl("^ Ighm*", results$Symbols),]
results <- results[!grepl("^ Ighv*", results$Symbols),]
results <- results[!grepl("^ Iglv*", results$Symbols),]
results <- results[!grepl("^ Igkc*", results$Symbols),]
results <- results[!grepl("^ Iglc*", results$Symbols),]


# Remove the whitespace from the symbols
results$Symbols <- trimws(results$Symbols)

# Save the results as a csv file
write.csv(results, file = "Outputs/002_KOSm_vs_WTSm_Outputs/DESeq2/KOSm_vs_WTSm_DESeq2_Results_IGs_Removed.csv")

# Unhash this line if you are editing figures
results <- read.csv("Outputs/002_KOSm_vs_WTSm_Outputs/DESeq2/KOSm_vs_WTSm_DESeq2_Results_IGs_Removed.csv")

#################################################################################
#----- Heatmap analysis

# Read in the log2(tpm) values
logTPM <- read.csv("Data_Files/TPM_Values/Log2_Transformed_TPM_Values_Delta14_project.csv") %>%
  column_to_rownames(var = "X")

# Get gene symbols
logTPM$Symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", rownames(logTPM)))

# Use the results file to order results by log2FC
results <- results[order(results$log2FoldChange, decreasing = TRUE),]

# Take the top and bottom 20
totalRows <- nrow(results)
bottom20 <- totalRows - 19
keep <- c(1:20, bottom20:totalRows)
topBottom <- results[keep,]$Symbols

# Subset the TPMS based on top and bottom 20
logTPM <- logTPM[logTPM$Symbols %in% topBottom,]
rownames(logTPM) <- logTPM$Symbols
logTPM$Symbols <- NULL

# Filter TPMs to only be samples for this comparison
meta <- metadata[metadata$Group == "C" | metadata$Group == "A",]$SampleID
logTPM <- logTPM[,colnames(logTPM) %in% meta]

# Get the dataframe in the proper order
logTPM <- logTPM[topBottom,]
# Calculate Z-score of the logTPMs
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate Zscore on log10 TPM data
mat <- t(apply(logTPM, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
n <- 7
m <- 7
ref <- "WT"
treatment <- "KO"

# Assign sample group colors
sample_colors <- c("steelblue2", "firebrick2")
names(sample_colors) <- c(ref, treatment)

#Set heatmap splitting pattern
hmSplit <- rep(1:2, c(n, m))

#Define the number of slices in the heatmap
slices <- n+m

#Create a heatmap annotation
anno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("steelblue2", "firebrick2"), fontisze = 14, fontface = "bold"), 
                     labels = c(ref, treatment),
                     labels_gp = gpar(col = "black", fontsize = 16, fontface = 2)),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)

#Heatmap for scaled data
scaled <- Heatmap(mat,
                  show_column_names = FALSE,
                  name = "Feature Scale",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  top_annotation = anno,
                  column_split = hmSplit,
                  column_title = NULL,
                  row_title = NULL,
                  row_names_gp = gpar(fontsize = 14))
tiff("Outputs/002_KOSm_vs_WTSm_Outputs/Figures/KOSm_vs_WTSm_TopBottom20_Zscaled_LogTPM_Heatmap.tiff", width = 8, height = 10, units = "in", res = 300)
draw(scaled)
dev.off()






