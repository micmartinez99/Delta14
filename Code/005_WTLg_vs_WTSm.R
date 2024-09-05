# The purpose of this code is to compare the WT large relative to the WT small.
# We already ran the model in script 004, so we can just load in the DESeq2 rds file
# And extract the comparison of intrest from there. 

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(PCAtools)
library(magick)
library(ComplexHeatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)

# Initialize function to generate output directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Create new output directory
opDir <- newDir("Outputs/005_WTLg_vs_WTSm_Outputs")

# Custom figures directory
figs <- newDir("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures")

# DESeq2 directory
ddsres <- newDir("Outputs/005_WTLg_vs_WTSm_Outputs/DESeq2")

# Rds files
rdsFiles <- newDir("Outputs/005_WTLg_vs_WTSm_Outputs/Rds_Files")

# GSEA outputs
gseares <- newDir("Outputs/005_WTLg_vs_WTSm_Outputs/GSEA")

################################################################################

#----- Exploratory data analysis (PCA)

# Read in the metadata
meta <- read.csv("Data_Files/Master_D14_Metadata.csv") %>%
  column_to_rownames(var = "SampleID")

# Read in the rds file
dds <- readRDS("Outputs/004_All_Group_Analysis/RDS_Files/All_Groups/AllGroups_Included_DESeq2.Rds")

# Get results names
resultsNames(dds)

# Extract the contrast to look at WT large vs WT small
# Here, we only need to look at Size Large vs small because our base level for 
# genotype is WT
results <- as.data.frame(results(dds, name = "Tumor_Size_Large_vs_Small"))

# Extract the normalized counts
normCounts <- as.data.frame(counts(dds, normalized = TRUE))

# Append results and counts
results <- cbind(results, normCounts)

# Extact the symbols
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

# Trim white space on the symbols
results$Symbols <- trimws(results$Symbols)

# Obtain Entrez IDs
results$Entrez <- mapIds(org.Mm.eg.db, 
                         keys = results$Symbols,
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

# Save the results as a csv
write.csv(results, file = "Outputs/005_WTLg_vs_WTSm_Outputs/DESeq2/WTLg_vs_WTSm_DESeq2_Results.csv")

# Filter the results for DEGs
res.filt <- results[results$baseMean > 50, ] 
res.filt <- res.filt[res.filt$padj <= 0.05 & abs(res.filt$log2FoldChange) > 1,]

# Save the filtered DEGs as a csv
write.csv(res.filt, file = "Outputs/005_WTLg_vs_WTSm_Outputs/DESeq2/WTLg_vs_WTSm_DEGs.csv")

################################################################################

# Variance stabilize transform the dds object
vsd <- vst(dds)

# Get gene symbols
symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", rownames(dds)))

# Define the patterns to remove
remove <- c("^Gm\\d+$", "Rik$", "Rik\\d+$", "^LOC", "AK\\d+$", "AY\\d+$",
            "^Igkv", "^Ighg", "^Igha", "^Ighd", "^Ighe", "^Ighm", "^Ighv", 
            "^Iglv", "^Igkc", "^Iglc")

# Create a logical vector to identify rows to keep
keep_genes <- !sapply(remove, function(pattern) grepl(pattern, symbols))

# Only keep the genes that do not match any of the patterns
keep_genes <- apply(keep_genes, 1, all)

# Subset the vsd object to keep only desired genes
vsd2 <- vsd
vsd2 <- vsd[keep_genes,]

# Sort the genes by decreasing variance
vsdassay <- sort(rowVars(assay(vsd2)), decreasing = TRUE)[1:500]

# Extract the names of the top 500 most variable features
varGenes <- names(vsdassay)

# Extract these genes from the vsd assay
pcaGenes <- assay(vsd2)
pcaGenes <- pcaGenes[rownames(pcaGenes) %in% varGenes,]
rownames(pcaGenes) <- gsub("^[^-]+-(.*)$", "\\1", rownames(pcaGenes))

# Keep only WT samples
meta$SampleID <- rownames(meta)
keepSamples <- meta[meta$Group == "B" | meta$Group == "A",]$SampleID
vsd3 <- vsd2[,colnames(vsd2) %in% keepSamples]

# Sort the genes by decreasing variance
vsdassay <- sort(rowVars(assay(vsd3)), decreasing = TRUE)[1:500]

# Extract the names of the top 500 most variable features
varGenes <- names(vsdassay)

# Extract these genes from the vsd assay
pcaGenes <- assay(vsd3)
pcaGenes <- pcaGenes[rownames(pcaGenes) %in% varGenes,]
rownames(pcaGenes) <- trimws(gsub("^[^-]+-(.*)$", "\\1", rownames(pcaGenes)))

# Run PCAtools
p2 <- pca(pcaGenes, metadata = colData(vsd3), removeVar = 0.1)

# Save PCA object as an Rds
saveRDS(p2, file = "Outputs/005_WTLg_vs_WTSm_Outputs/Rds_Files/WTLg_vs_WTSm_PCATools_Object.Rds")

# Create a subdirectory to hold the PCA analysis
pcaDir <- newDir("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures/PCA")

# Assess a biplot of just the genotype
class_biplot <- biplot(p2,
                       pointSize = 5,
                       ntopLoadings = 3,
                       showLoadings = TRUE,
                       sizeLoadingsNames = 5,
                       colby = "Class",
                       legendPosition = "bottom")
ggsave("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures/PCA/WTLg_vs_WTSm_Biplot.png", class_biplot, width = 8, height = 8)

################################################################################

#----- Heatmap analysis of DEGs
# Read in the log2(tpm) values
logTPM <- read.csv("Data_Files/TPM_Values/Log2_Transformed_TPM_Values_Delta14_project.csv") %>%
  column_to_rownames(var = "X")

# Get gene symbols
logTPM$Symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", rownames(logTPM)))

# Use the results file to order results by log2FC
res.filt <- res.filt[order(res.filt$log2FoldChange, decreasing = TRUE),]

# Get a vector of significant gene names
genesKeep <- res.filt$Symbols

# Subset the TPM file to only include these gene names
logTPM <- logTPM[logTPM$Symbols %in% genesKeep,]
rownames(logTPM) <- logTPM$Symbols
logTPM$Symbols <- NULL

# Filter TPMs to only be samples for this comparison
groups <- meta[meta$Group == "C" | meta$Group == "D",]$SampleID
logTPM <- logTPM[,colnames(logTPM) %in% groups]

# Get the dataframe in the proper order
logTPM <- logTPM[genesKeep,]

# Calculate Z-score of the logTPMs
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate Zscore on log10 TPM data
mat <- t(apply(logTPM, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
n <- 7
m <- 7
ref <- "WT Small"
treatment <- "WT Large"

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
                  name = "Z-Score",
                  cluster_rows = TRUE,
                  cluster_columns = FALSE,
                  top_annotation = anno,
                  column_split = hmSplit,
                  column_title = NULL,
                  row_title = NULL,
                  row_names_gp = gpar(fontsize = 14))
tiff("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures/WTLg_vs_WTSm_DEG_Heatmap.tiff", width = 8, height = 10, units = "in", res = 300)
draw(scaled)
dev.off()

################################################################################

# Order the genes and remove NA values
results <- results[order(results$log2FoldChange, decreasing = TRUE),]
results <- na.omit(results)

# Get a vector of fold changes and name them by their Entrez ID
genes <- results$log2FoldChange
names(genes) <- results$Entrez

# Remove any duplicate Entrez IDs and for sanity, order in decreasing order again
unique_entrez_genes <- names(genes[!duplicated(names(genes))])
unique_genes <- genes[unique_entrez_genes]
unique_genes <- sort(unique_genes, decreasing = TRUE)

# Check how many genes there are
length(unique_genes)

# Set seed for reproducibility
set.seed(03061999)

# Run GSEA for GO terms
GO <- gseGO(unique_genes,
            ont = "all",
            OrgDb = "org.Mm.eg.db",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            eps = 1e-300,
            verbose = TRUE,
            by = "fgsea",
            seed = TRUE)

# Simplify the  redundant GO terms
GO.simp <- clusterProfiler::simplify(GO, cutoff = 0.6, by = "p.adjust", select_fun = min)

# Save Rds
saveRDS(GO.simp, file = "Outputs/005_WTLg_vs_WTSm_Outputs/Rds_Files/WTLg_vs_WTSm_GSEA_GO.Rds")

# Save results as a dataframe and write to csv
GO.df <- as.data.frame(setReadable(GO.simp, "org.Mm.eg.db", "ENTREZID"))
write.csv(GO.df, file = "Outputs/005_WTLg_vs_WTSm_Outputs/GSEA/WTLarge_vs_WTSmall_GO_GSEA_Results.csv")

#----- GSEA visualizations

# Create output directory for GSEA figures
gseaFigs <- newDir("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures/GSEA_Figs")

# Calculate pairwise terms
pwt <- pairwise_termsim(GO.simp)
p1 <- treeplot(pwt)
ggsave("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures/GSEA_Figs/WTLarge_vs_WTSmall_GSEA_Treeplot.png", p1, width = 20, height = 10)


# Create a directory to hold barcode plots
newDir("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures/GSEA_Figs/BarCode_Plots")

# Order the gseGO dataframe by increasing enrichment score
GO.df <- GO.df[order(GO.df$enrichmentScore, decreasing = FALSE),]

# Get a vector of GO IDs
GOs <- unique(GO.df$ID)

# Itarate through the GO terms and plot
for (i in GOs) {
  term <- GO.df[GO.df$ID == i,]$Description
  print(term)
  
  # Format the term for fileNames
  termFormatted <- gsub(" ", "_", term)
  
  # Plot
  plot <- gseaplot2(GO.simp, geneSetID = i, pvalue_table = TRUE)
  
  # Set fileName
  fileName <- paste(termFormatted, ".png", sep = "")
  
  # Save
  ggsave(paste("Outputs/005_WTLg_vs_WTSm_Outputs/Custom_Figures/GSEA_Figs/BarCode_Plots", fileName, sep = "/"), plot, width = 10, height = 10)
}















