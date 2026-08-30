setwd("//wsl.localhost/Ubuntu-22.04/home/mda/RNAseq_project.nf/results.nf/counts")

counts_raw <- read.table("gene_counts.txt", header=TRUE, row.names=1, comment.char="#", check.names=FALSE)
counts_data <- counts_raw[, 6:ncol(counts_raw), drop=FALSE]
colnames(counts_data) <- gsub(".*(SRR[0-9]+).*", "\\1", colnames(counts_data))
head(counts_data)


n_samples <- ncol(counts_data)
print(paste("Total samples found:", n_samples))

counts_data <- counts_data[, order(colnames(counts_data))]
colnames(counts_data)

samples <- colnames(counts_data)
condition <- factor(c(rep("Control", 3), rep("Treated", 3)), levels = c("Control", "Treated"))
colData <- data.frame(condition = condition, row.names = samples)
print("=== colData Table ===")
print(colData)

write.csv(colData, "metadata.csv")
colData <- read.csv("metadata.csv", row.names = 1)
colData$condition <- factor(colData$condition, levels = c("Control", "Treated"))


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = colData, 
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "Treated", "Control"))
res_df <- as.data.frame(res)
summary(res)


library(org.Hs.eg.db)
res_df <- as.data.frame(res)
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

write.csv(res_df, "DESeq2_full_results.csv")

sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1.5)

print("=== TOP 10 BIOMARKERS ===")
head(sig_genes[order(sig_genes$padj), c("symbol", "log2FoldChange", "padj")], 10)


library(ggplot2)
res_df$expression[res_df$log2FoldChange > 0.585 & res_df$padj < 0.05] <- "UP"
res_df$expression[res_df$log2FoldChange < -0.585 & res_df$padj < 0.05] <- "DOWN"

#Volacano plot
ggplot(data = res_df[!is.na(res_df$padj), ], aes(x = log2FoldChange, y = -log10(padj), col = expression)) +
  geom_point(alpha = 0.6, size = 1.8) +
  theme_minimal() +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.585, 0.585), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  labs(title = "Differential Gene Expression (Volcano Plot)",
       subtitle = "Type 2 Diabtes vs Control",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value")


library(pheatmap)
vsd <- vst(dds, blind = FALSE)
top25_ids <- head(rownames(sig_genes[order(sig_genes$padj), ]), 25)

mat <- assay(vsd)[top25_ids, ]

symbols <- res_df[top25_ids, "symbol"]
rownames(mat) <- ifelse(is.na(symbols), top25_ids, symbols)

pheatmap(mat,
         scale = "row",
         clustering_distance_rows = "correlation",
         annotation_col = colData,
         main = "Top 25 Differentially Expressed Genes (Heatmap)")

library(clusterProfiler)
library(org.Hs.eg.db)
all_sig_genes <- rownames(res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ])
entrez_ids <- bitr(all_sig_genes,
                   fromType = "ENSEMBL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

ego <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.20)

ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
print(paste("Total Enriched Pathways Found:", nrow(as.data.frame(ego))))

# 5. Dotplot Render Karein
dotplot(ego, showCategory = 10, title = "Top Enriched Biological Processes")



# SSL error fix karne ke liye:
httr::set_config(httr::config(ssl_verifypeer = FALSE))

# Ab KEGG run karein
kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
dotplot(kegg, showCategory = 10, title = "Top Enriched KEGG Pathways")

ggsave("GO_Biological_Processes.png", width = 8, height = 6, dpi = 300)
ggsave("KEGG_Pathways.png", width = 8, height = 6, dpi = 300)
