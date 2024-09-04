# The purpose of this code is to examine the interaction effect in DESeq2.
# The major question here is why are small tumors in the KO being kept small?
# To examine this question, we can extract the interaction term when we specify a 
# more complicated model.

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
library(DESeq2)
library(PCAtools)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(magick)
library(ComplexHeatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)

# Initialize a new function to generate output directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Generate an output directory
opDir <- newDir("Outputs/004_All_Group_Analysis")

# Custom Figures director
figs <- newDir("Outputs/004_All_Group_Analysis/Custom_Figures")

################################################################################

# Read in the full non-normalized counts matrix
counts <- read.csv("Data_Files/Master_D14_Counts.csv") %>%
  column_to_rownames(var = "X")

# Read in the full metadata
meta <- read.csv("Data_Files/Master_D14_Metadata.csv") %>%
  column_to_rownames(var = "X")

# Make samples in the same order between counts and meta
sampleOrder <- rownames(meta)
counts <- counts[,sampleOrder]

# Sanity check that samples are not scrambled between counts and meta
all(colnames(counts) %in% rownames(meta))
all(colnames(counts) == rownames(meta))

# Prepare the DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Tumor_Size + Genotype + Tumor_Size:Genotype)

# Set a reference levels
dds$Genotype <- relevel(dds$Genotype, ref = "WT")
dds$Tumor_Size <- relevel(dds$Tumor_Size, ref = "Small")

# Run DESeq2
dds <- DESeq(dds)

# Check results
resultsNames(dds)

# Extract a snapshot of the results when comparing KOLarge vs KOSmall

# Explanation of contrast term: This contrast involves the main effect of tumor size (Large vs Small)
# and the interaction term captures the additional effect in the KO group.
# A + log2FC indicates higher gene expression in KOLg tumors
# A - log2FC indicates higher gene expression in KOSm tumors
snapshot <- results(dds, contrast = list(c("Tumor_Size_Large_vs_Small", "Tumor_SizeLarge.GenotypeKO"))) 
summary(snapshot, alpha = 0.05)

# Generate output directory for rds files
rds <- newDir("Outputs/004_All_Group_Analysis/RDS_Files")

# Extract the results
results <- as.data.frame(snapshot)
resultsFilt <- results[results$padj <= 0.05 & abs(results$log2FoldChange) > 1,]
resultsFilt <- na.omit(resultsFilt)

# Save dds object
saveRDS(dds, file = "AllGroups_Included_DESeq2.Rds")


################################################################################

#----- Prinicipal component analysis

# Variance stabilize transform the dds object
vsd <- vst(dds)

# Function to plot variances to determine a suitable number of features
determineVarFeatures <- function(vsd) {
  
  # Calculate gene expression level variance between samples
  var <- sort(rowVars(assay(vsd)), decreasing = TRUE)
  
  # Reset plotting window to 1 row vs 1 column
  par(mfrow = c(1, 1))
  
  # Plot the variance for genes across samples
  varPlot <- plot(var,
                  las = 1,
                  main = "Sample gene expression variance",
                  xlab = "Gene",
                  ylab = "Variance")
  
  # Add vertical lines as specific gene number indices
  abline(v = 1000, col = "red")
  abline(v = 500, col = "green")
  abline(v = 250, col = "blue")
  
  # Return values
  return(varPlot)
}

# Run determineVarFeatures function
determineVarFeatures(vsd)

# Run PCA and return data for custom plotting
PCA <- plotPCA(vsd, intgroup = c("Class"), returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))

# Set custom colors
custom_colors <- c("KO_Large" = "cornflowerblue",
                   "KO_Small" = "firebrick",
                   "WT_Large" = "cyan4",
                   "WT_Small" = "goldenrod")

# Plot
ggplot(PCA, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 4) + 
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(x = "PC1: 43%",
       y = "PC2: 14%",
       color = "",
       fill = "") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = "bottom")


#----- PCA analysis, extended
pcaDir <- newDir("Outputs/004_All_Group_Analysis/Custom_Figures/PCA_Analysis")
byGeno <- newDir("Outputs/004_All_Group_Analysis/Custom_Figures/PCA_Analysis/All_Groups_By_Genotype")
PCArds <- newDir("Outputs/004_All_Group_Analysis/RDS_Files")
PCArdsSub <- newDir("Outputs/004_All_Group_Analysis/RDS_Files/All_Groups")
bySize <- newDir("Outputs/004_All_Group_Analysis/Custom_Figures/PCA_Analysis/All_Groups_By_Size")
byClass <- newDir("Outputs/004_All_Group_Analysis/Custom_Figures/PCA_Analysis/All_Groups_By_Group")

# Remove IGs and unknown genes
# Remove unknown genes 
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

# Run PCAtools
p <- pca(pcaGenes, metadata = colData(dds), removeVar = 0.1)

# Save the PCAtools object
saveRDS(p, file = "Outputs/004_All_Group_Analysis/RDS_Files/All_Groups/PCATools_IGs_Removed.rds")

# Assess a biplot of just the genotype
genotype_biplot <- biplot(p,
                          pointSize = 5,
                          ntopLoadings = 3,
                          showLoadings = TRUE,
                          sizeLoadingsNames = 7,
                          colby = "Genotype",
                          legendPosition = "bottom")
ggsave("Outputs/004_All_Group_Analysis/Custom_Figures/PCA_Analysis/All_Groups_By_Genotype/AllGroups_IGs_Removed_By_Genotype_Biplot.png", genotype_biplot, width = 8, height = 8)

# Assess biplot based on group
group_biplot <- biplot(p,
                       pointSize = 5,
                       ntopLoadings = 3,
                       showLoadings = TRUE,
                       sizeLoadingsNames = 7,
                       colby = "Class",
                       legendPosition = "bottom")
ggsave("Outputs/004_All_Group_Analysis/Custom_Figures/PCA_Analysis/All_Groups_By_Group/AllGroups_IGs_Removed_By_Group_Biplot.png", group_biplot, width = 8, height = 8)

keepSamples <- meta[meta$Group == "C" | meta$Group == "D",]$SampleID
vsd3 <- vsd2[,colnames(vsd2) %in% keepSamples]

# Sort the genes by decreasing variance
vsdassay <- sort(rowVars(assay(vsd3)), decreasing = TRUE)[1:500]

# Extract the names of the top 500 most variable features
varGenes <- names(vsdassay)

# Extract these genes from the vsd assay
pcaGenes <- assay(vsd3)
pcaGenes <- pcaGenes[rownames(pcaGenes) %in% varGenes,]
rownames(pcaGenes) <- gsub("^[^-]+-(.*)$", "\\1", rownames(pcaGenes))

# Run PCAtools
p2 <- pca(pcaGenes, metadata = colData(vsd3), removeVar = 0.1)

# Save the PCAtools object
saveRDS(p2, file = "Outputs/004_All_Group_Analysis/RDS_Files/All_Groups/PCATools_IGs_Removed_KO_vs_KO.rds")

# Assess a biplot of just the genotype
class_biplot <- biplot(p2,
                          pointSize = 5,
                          ntopLoadings = 3,
                          showLoadings = TRUE,
                          sizeLoadingsNames = 5,
                          colby = "Class",
                          legendPosition = "bottom")
ggsave("Outputs/004_All_Group_Analysis/Custom_Figures/PCA_Analysis/KO_vs_KO_Biplot.png", class_biplot, width = 8, height = 8)

################################################################################

# Initialize a new dir for DESeq2 results
newDir("Outputs/004_All_Group_Analysis/DESeq2_Results")

# Extract the results for KO Large vs KO Small
results <- as.data.frame(results(dds, contrast = list(c("Tumor_Size_Large_vs_Small", "Tumor_SizeLarge.GenotypeKO"))))

# Extract the normalized counts
normCounts <- as.data.frame(counts(dds, normalized = TRUE))

# Combine the results
results <- cbind(results, normCounts)

# Omit genes with no padj value
results <- na.omit(results)

# Filter out IG genes and other unknowns
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

# Trim whitespace from symbols
results$Symbols <- trimws(results$Symbols)

# Obtain Entrez IDs
results$Entrez <- mapIds(org.Mm.eg.db, 
                         keys = results$Symbols,
                         column = "ENTREZID",
                         keytype = "SYMBOL")

# Save as a csv
write.csv(results, file = "Outputs/004_All_Group_Analysis/DESeq2_Results/KOLarge_vs_KOSmall_DESeq2_Results.csv")

# Filter the results for significant genes
results <- results[results$baseMean > 50,]
res.filt <- results[results$padj <= 0.05 & abs(results$log2FoldChange) > 1,]

# Save the significant DEGs as a csv
write.csv(res.filt, file = "Outputs/004_All_Group_Analysis/DESeq2_Results/KOLarge_vs_KOSmall_DEGs.csv")

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
ref <- "KO Small"
treatment <- "KO Large"

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
tiff("Outputs/004_All_Group_Analysis/Custom_Figures/KOLarge_vs_KOSmall_DEG_Heatmap.tiff", width = 8, height = 10, units = "in", res = 300)
draw(scaled)
dev.off()

################################################################################

#----- GSEA

# Generate a GSEA output dir
newDir("Outputs/004_All_Group_Analysis/GSEA_Results")

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

# Save gseaResult object as a .Rds object
saveRDS(GO.simp, file = "Outputs/004_All_Group_Analysis/RDS_Files/KOLarge_vs_KOSmall_GO_GSEA.rds")

# Save results as a dataframe and write to csv
GO.df <- as.data.frame(setReadable(GO.simp, "org.Mm.eg.db", "ENTREZID"))
write.csv(GO.df, file = "Outputs/004_All_Group_Analysis/GSEA_Results/KOLarge_vs_KOSmall_GO_GSEA_Results.csv")


#----- GSEA visualizations

# Create output directory for GSEA figures
gseaFigs <- newDir("Outputs/004_All_Group_Analysis/Custom_Figures/GSEA_Figs")

# Calculate pairwise terms
pwt <- pairwise_termsim(GO.simp)
p1 <- treeplot(pwt)
ggsave("Outputs/004_All_Group_Analysis/Custom_Figures/GSEA_Figs/KOLarge_vs_KOSmall_GSEA_Treeplot.png", p1, width = 20, height = 10)

# Upset plot for GSEA
upsetplot(GO.simp) 

# Create a directory to hold barcode plots
newDir("Outputs/004_All_Group_Analysis/Custom_Figures/GSEA_Figs/BarCode_Plots")

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
  ggsave(paste("Outputs/004_All_Group_Analysis/Custom_Figures/GSEA_Figs/BarCode_Plots", fileName, sep = "/"), plot, width = 10, height = 10)
}


#----- Let's see what genes are involved in some of the interesting GO terms
# Get a vector of fold changes and name them by their Entrez ID
genes <- results$log2FoldChange
names(genes) <- results$Symbols

# Remove any duplicate Symbols and for sanity, order in decreasing order again
unique_entrez_genes <- names(genes[!duplicated(names(genes))])
unique_genes <- genes[unique_entrez_genes]
unique_genes <- sort(unique_genes, decreasing = TRUE)


# Make GO results readable
GO.readable <- setReadable(GO.simp, "org.Mm.eg.db", "ENTREZID")

# Plot CNET
CNETPlot <- cnetplot(GO.readable, 
         foldChange = unique_genes,
         circular = TRUE, 
         colorEdge = TRUE,
         node_label = "gene",
         showCategory = c("carboxylesterase activity",
                          "positive regulation of leukocyte tethering or rolling",
                          "TAP binding",
                          "positive regulation of CD8-positive, alpha-beta T cell proliferation",
                          "innate immune response in mucosa",
                          "antimicrobial humoral immune response mediated by antimicrobial peptide"))
ggsave("Outputs/004_All_Group_Analysis/Custom_Figures/GSEA_Figs/KOLarge_vs_KOSmall_CNET.png", CNETPlot, width = 10, height = 10)


