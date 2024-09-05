# Delta14

# Analysis Methods
-000_TPM_Generation.R
The purpose of this script is to take the master expression matrix (all samples and all mapped reads) and convert these counts to transcripts per million (TPM) values. To successfully reproduce this code, the file `Data_Files/Mmu_Gene_Lengths/Mus_musculus_exon_lengths_in_BP.csv" is provided. As the file name implies, these gene lengths are in terms of basepairs and not kilo-basepairs. This code handles the conversion from bp to kbp, then converts the counts matrix to TPM values. This script will output 2 files: Non-transformed TPM values and log2 transformed TPM values (for heatmap visualizations.)

-001_KOLg_WTLg.R
This code handles the inter-genotype comparison for the large tumor samples. For this comparison, we are only interested in directly comparing these two groups so the design is specified as a grouping variable. After running DESeq2, immunoglobulins and unknown genes were removed since this particular KO model directly influences IG gene production, causing immunoglobulins to act as "noise" in our DEG results. 

-002_KOSm_vs_WTSm.R
Same logic as in script 001, but this comparison is the inter-genotype comparison for small tumors. 

-003_VennDiagram
This code was used to generate the number of differentially expressed genes in each comparison. (I ended up needing to go back to this code after running script 004 to include the KO vs KO comparison to generate an upset plot.)

-004_KOLg_vs_KOSm.R
This code deals with the within-genotype comparison for the KO model. Here, we are interested in assessing the question, what gene expression changes could be keeping the small tumors from progressing to large tumors? As above, after running DESeq2, immunoglobulins and unknown genes were removed. We probe this question by using principal component analysis for all groups first, then only for the KOLg and KOSm groups. After filtering for genes with baseMean > 50, padj < 0.05, and abs(log2FC) > 1, we identified a small subset of genes differentially expressed between the two groups. 
Gene set enrichment analysis was used to assess coordinated changes in the list of genes output by DESeq2.




