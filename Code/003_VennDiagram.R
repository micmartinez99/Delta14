# The purpose of this code is to analyze changes in differentially expressed genes
# The broad question is what is keeping tumors small in the KO model?
# To assess this, we ran two comparisons: KO vs WT for large and KO vs WT for small
# We want to visualize what is UP in the small KO that is also DOWN in the large KO

# Clear the environment
rm(list = ls())

# Load in libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggvenn)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Create output directories
opDir <- "Outputs/003_VennDiagram_Outputs"
if (!dir.exists(opDir)) {
  dir.create(opDir)
}

# Create a subdirectory for figures
figDir <- paste(opDir, "Figures", sep = "/")
if (!dir.exists(figDir)) {
  dir.create(figDir)
}

#################################################################################

# Read in the results for the large comparisons
large <- read.csv("Outputs/001_KOLg_vs_WTLg_Outputs/DESeq2/KOLg_vs_WTLg_DESeq2_Results_IGs_Removed.csv")
large <- large %>%
  column_to_rownames(var = "X")

# Read in the results for the small comparison
small <- read.csv("Outputs/002_KOSm_vs_WTSm_Outputs/DESeq2/KOSm_vs_WTSm_DESeq2_Results_IGs_Removed.csv")
small <- small %>%
  column_to_rownames(var = "X")

# Get only the genes in common between the two comparisons
common <- intersect(rownames(large), rownames(small))

# Filter both data-frames using the "common" vector
large <- large[rownames(large) %in% common,]
small <- small[rownames(small) %in% common,]

# Genes that are UPREGULATED in the small
smallUp <- small[small$log2FoldChange > 1 & small$padj < 0.05, ]

# Genes that are UPREGULATED in the large KO
largeUp <- large[large$log2FoldChange > 1 & large$padj < 0.05, ]

# How many of these genes are shared between smallUp and largeDown
intersect(rownames(smallUp), rownames(largeUp))

# Make a venn diagram
a <- list(smUp = smallUp$Symbols,
          lgUp = largeUp$Symbols)

# DOWN in treatment groups
A <- ggvenn(a, c("smUp", "lgUp"), show_percentage = FALSE,
               fill_color = c("red", "cyan"), text_size = 12, set_name_size = 10, auto_scale = TRUE) +
  theme(text = element_text(family = "Times New Roman"))
ggsave("Outputs/003_VennDiagram_Outputs/Figures/KO_Small_KO_Large_Upregulation_Venn_Diagram.tiff", A, width = 10, height = 10, dpi = 300)

################################################################################

# Now let's look at the genes that are DOWN in KO small and DOWN in KO large
smallDn <- small[small$log2FoldChange < -1, ]

# Genes that are DOWNREGULATED in the large KO
largeDn <- large[large$log2FoldChange < -1, ]

# Make a venn diagram
b <- list(smDn = smallDn$Symbols,
          lgDn = largeDn$Symbols)

# DOWN in treatment groups
B <- ggvenn(b, c("smDn", "lgDn"), show_percentage = FALSE,
            fill_color = c("red", "cyan"), text_size = 12, set_name_size = 10, auto_scale = TRUE) +
  theme(text = element_text(family = "Times New Roman"))


# Now let's look at the genes that are UP in KO small and DOWN in KO large

# Make a venn diagram
c <- list(smUp = smallUp$Symbols,
          lgDn = largeDn$Symbols)

# DOWN in treatment groups
C <- ggvenn(c, c("smUp", "lgDn"), show_percentage = FALSE,
            fill_color = c("red", "cyan"), text_size = 12, set_name_size = 10, auto_scale = TRUE) +
  theme(text = element_text(family = "Times New Roman"))


################################################################################
# Extract the list of genes that are UP in KO small
smUp <- a[["smUp"]]

# Remove the genes that are shared between smUp and lgUP
same <- intersect(smUp, largeUp$Symbols)
smUp <- smallUp[smallUp$Symbols %in% smUp,]
smUp <- smUp[!(smUp$Symbols) %in% same, ]

# Get gene Information
smUp$Desc <- mapIds(org.Mm.eg.db, key = smUp$Symbols,
                         column = "GENENAME", keytype = "SYMBOL",
                         multiVals = "first")
# Set the rownames
rownames(smUp) <- smUp$Symbols
smUp$Symbols <- NULL
smUp <- smUp[,c(1,2,6,21)]
smUp$Symbols <- rownames(smUp)
smUp$Group <- c("Up regulated in KO Small")

# Get the genes that are shared
smallSame <- small[small$Symbols %in% same,]

# Get gene Information
smallSame$Desc <- mapIds(org.Mm.eg.db, key = smallSame$Symbols,
                    column = "GENENAME", keytype = "SYMBOL",
                    multiVals = "first")
rownames(smallSame) <- smallSame$Symbols
smallSame$Symbols <- NULL
smallSame <- smallSame[,c(1,2,6,21)]
smallSame$Symbols <- rownames(smallSame)
smallSame$Group <- "Up regulated in KO Small and Large"

# Get only the genes exclusively upregulated in the large KO
lgUp <- a[["lgUp"]]
lgUp <- largeUp[largeUp$Symbols %in% lgUp,]
lgUp <- lgUp[!(lgUp$Symbols) %in% same, ]

# Get gene Information
lgUp$Desc <- mapIds(org.Mm.eg.db, key = lgUp$Symbols,
                    column = "GENENAME", keytype = "SYMBOL",
                    multiVals = "first")

# Set the rownames
rownames(lgUp) <- lgUp$Symbols
lgUp$Symbols <- NULL
lgUp <- lgUp[,c(1,2,6,21)]
lgUp$Symbols <- rownames(lgUp)
lgUp$Group <- c("Upregulated in KO Large")

# Combine the results 
results <- rbind(rbind(smUp, smallSame), lgUp)
write.csv(results, file = "Outputs/003_VennDiagram_Outputs/UpInSmallKO_UpInLargeKO_Venn_Results.csv")

################################################################################

# Upset plot showing the common genes
library(UpSetR)

sameDown <- intersect(largeDn$Symbols, smallDn$Symbols)

# Read in the KO vs KO results
KOvKO <- read.csv("Outputs/004_All_Group_Analysis/DESeq2_Results/KOLarge_vs_KOSmall_DEGs.csv")
KOvKOsmall <- KOvKO[KOvKO$log2FoldChange < -1,]
KOvKOlarge <- KOvKO[KOvKO$log2FoldChange > 1,]


# Upset plot expression matrix
expressionInput <- c(KO.Large = 471, # Up in KO Large (KO large vs WT large)
                     KO.Small = 340, # Up in KO Small (KO small vs WT small)
                     WT.Large = 296, # Down in KO Large (KO large vs WT large)
                     WT.Small = 494, # Down in KO Small (KO Small vs WT Small)
                     KO.Small_Specific = 16, # Genes down regulated in KOLg vs KOSmall,
                     KO.Large_Specific = 10, # Genes up regulated in KOLg vs KOSmall
                     `KO.Large&KO.Small` = 235,
                     `WT.Large&WT.Small` = 238)


# Plot Upset plot
upsetPlot <- upset(fromExpression(expressionInput),
                   order.by = "freq",
                   decreasing = TRUE)

tiff("Outputs/004_All_Group_Analysis/Custom_Figures/UpsetPlot.tiff", width = 8, height = 10, units = "in", res = 300)
draw(upsetPlot)
dev.off()














