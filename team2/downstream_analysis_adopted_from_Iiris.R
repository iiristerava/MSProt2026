library(limma)
library(pheatmap)
library(tidyverse)
library(knitr)
library(ggrepel)
library(arrow)
library(data.table)
library(pcaMethods)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(grid)
library(ggVennDiagram)

# Output Directory
out_dir <- "/Users/mikhaila/Desktop/MSProt2026/team2/plots_1506"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

parquet_file <- "/Users/mikhaila/Desktop/MSProt2026/team2/transfer_421909_files_47374ab3/report.parquet"
fasta_file <- "/Users/mikhaila/Desktop/MSProt2026/team2/transfer_421909_files_47374ab3/contaminants.fasta"


data <- open_dataset(parquet_file)

# Contaminants
cont_fasta <- readLines(fasta_file)
cont_headers <- cont_fasta[grepl("^>", cont_fasta)] 
cont_ids <- sub("^>([^ ]+).*", "\\1", cont_headers) 

# Proteotypic Precursors Only
data_precursor <- data %>%
  filter(Proteotypic == TRUE) %>%
  select(Precursor.Id, Run, Precursor.Normalised, Protein.Ids) %>%
  collect()

# Filter Singletons
tab <- table(data_precursor$Precursor.Id)
data_precursor <- data_precursor %>%
  group_by(Precursor.Id) %>%
  filter(n() > 1) %>%
  ungroup()

# Fast Contaminant Removal
search_pattern <- paste(cont_ids, collapse = "|")
is_contaminant <- str_detect(data_precursor$Protein.Ids, search_pattern)
is_contaminant[is.na(is_contaminant)] <- FALSE
data_precursor <- data_precursor[!is_contaminant, ]

# Reshape to Matrix
dt <- as.data.table(data_precursor)
matrix_data <- dcast(dt, Precursor.Id ~ Run, value.var = "Precursor.Normalised")

# Missing Value Filter (30% Rule)
matrix_data_filtered <- matrix_data[rowSums(!is.na(matrix_data)) >= 0.3 * ncol(matrix_data), ]
matrix_data_filtered <- as.data.frame(matrix_data_filtered)
rownames(matrix_data_filtered) <- matrix_data_filtered$Precursor.Id
matrix_data_filtered$Precursor.Id <- NULL

# Annotation Table
annotation <- data %>%
  filter(Proteotypic == TRUE) %>%
  select(Precursor.Id, Stripped.Sequence, Protein.Names, Genes) %>%
  collect() %>%
  unique()


# Log2 & Median Normalization
log_data <- log2(as.matrix(matrix_data_filtered) + 1)
sample_medians <- apply(log_data, 2, median, na.rm = TRUE)
normalized_data <- sweep(log_data, 2, sample_medians)

# Metadata Assignment
sample_names <- colnames(normalized_data)

group <- ifelse(grepl("CTRL", sample_names) & !grepl("cDDP", sample_names), "CTRL_untreated",
                ifelse(grepl("KO", sample_names) & !grepl("cDDP", sample_names), "KO_untreated",
                       ifelse(grepl("CTRL", sample_names) & grepl("cDDP.*06h", sample_names) & !grepl("JNK-IN-8", sample_names), "CTRL_cDDP_6h",
                              ifelse(grepl("KO", sample_names) & grepl("cDDP.*06h", sample_names), "KO_cDDP_6h",
                                     ifelse(grepl("CTRL", sample_names) & grepl("cDDP.*24h", sample_names) & !grepl("JNK-IN-8", sample_names), "CTRL_cDDP_24h",
                                            ifelse(grepl("KO", sample_names) & grepl("cDDP.*24h", sample_names), "KO_cDDP_24h",
                                                   ifelse(grepl("JNK-IN-8", sample_names) & grepl("06h", sample_names), "CTRL_JNK_cDDP_6h",
                                                          ifelse(grepl("JNK-IN-8", sample_names) & grepl("24h", sample_names), "CTRL_JNK_cDDP_24h", "Other"))))))))
group <- factor(group)

prep_group <- ifelse(grepl("20260331", sample_names), "Group0_20260331",
                     ifelse(grepl("Group1", sample_names), "Group1",
                            ifelse(grepl("Group2", sample_names), "Group2",
                                   ifelse(grepl("Group3", sample_names), "Group3",
                                          ifelse(grepl("Group4", sample_names), "Group4", "Other")))))
prep_group <- factor(prep_group)

# Batch Correction (Iiris's Method: Create modified matrix for ALL downstream steps)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

normalized_data_corrected <- removeBatchEffect(normalized_data, batch = prep_group, design = design)


pca_matrix <- normalized_data_corrected[complete.cases(normalized_data_corrected), ]
pca_res <- prcomp(t(pca_matrix), scale. = TRUE)

pca_data <- as.data.frame(pca_res$x)
pca_data$File.Name <- rownames(pca_data)

meta_data <- data.frame(File.Name = sample_names, Group = group)
pca_full <- merge(pca_data, meta_data, by="File.Name")
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

custom_colors <- c(
  "KO_cDDP_24h"       = "#56C1FF",   
  "KO_cDDP_6h"        = "#ED9A8D",   
  "KO_untreated"      = "#CF9BFA",   
  "CTRL_cDDP_24h"     = "#A5CA39",   
  "CTRL_cDDP_6h"      = "#F378E0",   
  "CTRL_JNK_cDDP_24h" = "#5CD9D7",   
  "CTRL_JNK_cDDP_6h"  = "#DEAE3B",   
  "CTRL_untreated"    = "#6CD184"    
)

p_custom <- ggplot(pca_full, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = custom_colors) + 
  theme_minimal() +
  labs(title = "PCA of Proteomic Profiles (Batch Corrected)",
       x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)"))

print(p_custom)
ggsave(filename = file.path(out_dir, "PCA_Plot_Combined.pdf"), plot = p_custom, width = 10, height = 7)

# Feeding the removeBatchEffect matrix directly into lmFit (Iiris's method)
fit <- lmFit(normalized_data_corrected, design)

contrasts_matrix <- makeContrasts(
  Mechanistic_Baseline      = "KO_untreated - CTRL_untreated",
  Mechanistic_Stress        = "KO_cDDP_24h - CTRL_cDDP_24h",
  Mechanistic_Stress_1      = "KO_cDDP_6h - CTRL_cDDP_6h", 
  Translational_Inhibitor   = "CTRL_JNK_cDDP_24h - CTRL_cDDP_24h",
  Translational_Inhibitor_1 = "CTRL_JNK_cDDP_6h - CTRL_cDDP_6h",
  Chemical_Genetic          = "CTRL_JNK_cDDP_24h - KO_cDDP_24h",
  Chemical_Genetic_1        = "CTRL_JNK_cDDP_6h - KO_cDDP_6h",
  Biomarker_Resistance      = "CTRL_cDDP_24h - CTRL_untreated",
  levels = design
)

fit_contrasts <- contrasts.fit(fit, contrasts_matrix)
fit_contrasts <- eBayes(fit_contrasts)

contrast_names <- colnames(fit_contrasts$contrasts)

for (contrast in contrast_names) {
  cat("\n\n### Contrast Group:", contrast, "\n")
  
  # --- 1. Prepare DEA Results Table ---
  res <- topTable(fit_contrasts, coef = contrast, number = Inf, adjust.method = "fdr")
  res$Precursor.Id <- rownames(res)
  
  # CRITICAL FIX: Sort before distinct to ensure lowest p-value label is kept
  res <- merge(res, annotation, by = "Precursor.Id", all.x = TRUE)
  res <- res[order(res$P.Value), ] 
  
  res$PlotLabel <- ifelse(is.na(res$Genes) | res$Genes == "", res$Precursor.Id, res$Genes)
  res_unique <- res %>% distinct(PlotLabel, .keep_all = TRUE)
  
  # --- 2. Volcano Plot ---
  to_label <- res_unique %>% filter(adj.P.Val < 0.01 & abs(logFC) > 1)
  if (nrow(to_label) > 0) {
    to_label <- to_label[order(to_label$adj.P.Val / abs(to_label$logFC)), ]
    top_15_genes <- head(to_label$PlotLabel, 15)
  } else {
    top_15_genes <- head(res_unique[order(res_unique$adj.P.Val), "PlotLabel"], 15)
  }
  
  v_plot <- EnhancedVolcano(res_unique, lab = res_unique$PlotLabel, selectLab = top_15_genes, 
                            x = 'logFC', y = 'adj.P.Val', pCutoff = 0.01, FCcutoff = 1.0, 
                            pointSize = 1.8, labSize = 3.0, drawConnectors = TRUE, widthConnectors = 0.4,
                            title = paste("Volcano:", str_replace_all(contrast, "_", " ")), legendPosition = 'bottom')
  print(v_plot)
  ggsave(filename = file.path(out_dir, paste0("Volcano_", contrast, ".pdf")), plot = v_plot, width = 9, height = 8)
  
  # --- 3. Annotated Heatmaps (Both Versions) ---
  
  # SAFETY FIX: Isolate only precursors that have NO missing values across the matrix
  complete_ids <- rownames(normalized_data_corrected)[complete.cases(normalized_data_corrected)]
  res_heatmap <- res_unique %>% filter(Precursor.Id %in% complete_ids)
  
  # Now safely select the top 50 features from the complete cases
  sig_features <- res_heatmap %>% filter(adj.P.Val < 0.05) %>% head(50)
  if(nrow(sig_features) < 15) { sig_features <- head(res_heatmap, 50) }
  
  if (nrow(sig_features) > 5) {
    involved_groups <- names(which(contrasts_matrix[, contrast] != 0))
    
    # Heatmap A: WITHOUT SUBSETTING (All samples shown)
    heatmap_matrix_all <- normalized_data_corrected[sig_features$Precursor.Id, , drop = FALSE]
    rownames(heatmap_matrix_all) <- sig_features$PlotLabel
    
    # Second safety check: remove any rows where variance is 0 (prevents scaling crashes)
    valid_rows_all <- apply(heatmap_matrix_all, 1, var) > 0
    heatmap_matrix_all <- heatmap_matrix_all[valid_rows_all, , drop = FALSE]
    
    if(nrow(heatmap_matrix_all) > 2) {
      annotation_col_all <- data.frame(Experimental_Condition = group, Comparison_Status = ifelse(group %in% involved_groups, "Compared", "Not Compared"))
      rownames(annotation_col_all) <- colnames(normalized_data_corrected)
      annotation_row <- data.frame(Expression = ifelse(sig_features$logFC[valid_rows_all] > 0, "Upregulated", "Downregulated"))
      rownames(annotation_row) <- rownames(heatmap_matrix_all)
      ann_colors_all = list(Experimental_Condition = custom_colors, Expression = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8"), Comparison_Status = c("Compared" = "black", "Not Compared" = "#E0E0E0"))
      
      ph_all <- pheatmap(heatmap_matrix_all, scale = "row", clustering_method = "ward.D2", annotation_col = annotation_col_all, annotation_row = annotation_row, annotation_colors = ann_colors_all, show_colnames = FALSE, fontsize_row = 7, silent = TRUE)
      
      pdf(file = file.path(out_dir, paste0("Heatmap_AllSamples_", contrast, ".pdf")), width = 10, height = 9)
      grid.newpage()
      grid.draw(ph_all$gtable)
      invisible(dev.off())
    }
    
    # Heatmap B: WITH SUBSETTING (Only compared samples shown)
    samples_to_keep <- rownames(annotation_col_all)[annotation_col_all$Comparison_Status == "Compared"]
    heatmap_matrix_sub <- normalized_data_corrected[sig_features$Precursor.Id, samples_to_keep, drop = FALSE]
    rownames(heatmap_matrix_sub) <- sig_features$PlotLabel
    
    # Recalculate variance for the specific subset to ensure no flat-lines
    valid_rows_sub <- apply(heatmap_matrix_sub, 1, var) > 0
    heatmap_matrix_sub <- heatmap_matrix_sub[valid_rows_sub, , drop = FALSE]
    
    if(nrow(heatmap_matrix_sub) > 2) {
      annotation_col_sub <- data.frame(Experimental_Condition = group[group %in% involved_groups])
      rownames(annotation_col_sub) <- samples_to_keep
      
      annotation_row_sub <- data.frame(Expression = ifelse(sig_features$logFC[valid_rows_sub] > 0, "Upregulated", "Downregulated"))
      rownames(annotation_row_sub) <- rownames(heatmap_matrix_sub)
      ann_colors_sub = list(Experimental_Condition = custom_colors, Expression = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8"))
      
      ph_sub <- pheatmap(heatmap_matrix_sub, scale = "row", clustering_method = "ward.D2", annotation_col = annotation_col_sub, annotation_row = annotation_row_sub, annotation_colors = ann_colors_sub, show_colnames = FALSE, fontsize_row = 7, silent = TRUE)
      
      pdf(file = file.path(out_dir, paste0("Heatmap_Subset_", contrast, ".pdf")), width = 10, height = 9)
      grid.newpage()
      grid.draw(ph_sub$gtable)
      invisible(dev.off())
    }
  }
}
for (contrast_name in contrast_names) {
  cat("\n### Individual Pathway Enrichment:", contrast_name, "\n")
  
  res <- topTable(fit_contrasts, coef = contrast_name, number = Inf, adjust.method = "fdr")
  res$Precursor.Id <- rownames(res)
  
  res <- merge(res, annotation, by = "Precursor.Id", all.x = TRUE)
  res <- res[order(res$P.Value), ]
  
  sig_res <- res %>% filter(adj.P.Val < 0.05 & !is.na(Genes) & Genes != "") %>% distinct(Genes, .keep_all = TRUE)
  
  if (nrow(sig_res) > 0) {
    conv <- bitr(sig_res$Genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (nrow(conv) > 0) {
      
      ego <- enrichGO(gene = conv$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        p_go <- dotplot(ego, showCategory = 10, title = paste("GO Enrichment:", str_replace_all(contrast_name, "_", " ")))
        pdf(file = file.path(out_dir, paste0("GO_", contrast_name, ".pdf")), width = 10, height = 8)
        print(p_go)
        invisible(dev.off())
      }
      
      ekegg <- enrichKEGG(gene = conv$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        p_kegg <- dotplot(ekegg, showCategory = 10, title = paste("KEGG Enrichment:", str_replace_all(contrast_name, "_", " ")))
        pdf(file = file.path(out_dir, paste0("KEGG_", contrast_name, ".pdf")), width = 10, height = 8)
        print(p_kegg)
        invisible(dev.off())
      }
    }
  }
}

get_sig_genes <- function(contrast_name) {
  res <- topTable(fit_contrasts, coef = contrast_name, number = Inf, adjust.method = "fdr")
  res$Precursor.Id <- rownames(res)
  
  res <- merge(res, annotation, by = "Precursor.Id", all.x = TRUE)
  res <- res[order(res$P.Value), ]
  
  sig_res <- res %>% filter(adj.P.Val < 0.05 & !is.na(Genes) & Genes != "")
  return(unique(sig_res$Genes))
}

list_6h <- list(`c-Jun KO (6h)` = get_sig_genes("Mechanistic_Stress_1"), `JNK-IN-8 (6h)` = get_sig_genes("Translational_Inhibitor_1"))
list_24h <- list(`c-Jun KO (24h)` = get_sig_genes("Mechanistic_Stress"), `JNK-IN-8 (24h)` = get_sig_genes("Translational_Inhibitor"))

venn_6h <- ggVennDiagram(list_6h, label_alpha = 0) + scale_fill_gradient(low = "#F4FAFE", high = "#56C1FF") + theme(legend.position = "none") + labs(title = "Gene Overlap at 6h")
print(venn_6h)
ggsave(filename = file.path(out_dir, "Venn_Overlap_6h.pdf"), plot = venn_6h, width = 7, height = 5)

venn_24h <- ggVennDiagram(list_24h, label_alpha = 0) + scale_fill_gradient(low = "#F4FAFE", high = "#A5CA39") + theme(legend.position = "none") + labs(title = "Gene Overlap at 24h")
print(venn_24h)
ggsave(filename = file.path(out_dir, "Venn_Overlap_24h.pdf"), plot = venn_24h, width = 7, height = 5)


# Convert Symbols directly to Entrez
convert_to_entrez <- function(gene_symbols) {
  if(length(gene_symbols) == 0) return(character(0))
  conv <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(conv$ENTREZID)
}

entrez_lists <- list(
  `KO 6h` = convert_to_entrez(list_6h[[1]]), `Chem 6h` = convert_to_entrez(list_6h[[2]]),
  `KO 24h` = convert_to_entrez(list_24h[[1]]), `Chem 24h` = convert_to_entrez(list_24h[[2]])
)
entrez_lists <- entrez_lists[lengths(entrez_lists) > 0]

# GO Biological Process
comp_GO <- compareCluster(geneCluster = entrez_lists, fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
if (!is.null(comp_GO) && nrow(as.data.frame(comp_GO)) > 0) {
  
  # CHANGED: showCategory = 5 (5 per group = max 20 total on the Y-axis)
  p_comp_GO <- dotplot(comp_GO, showCategory = 5, label_format = 45, title = "Pathway Intersection (GO)") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  pdf(file = file.path(out_dir, "Pathway_Intersection_GO.pdf"), width = 10, height = 8)
  print(p_comp_GO)
  invisible(dev.off())
  print(p_comp_GO)
}

# KEGG Pathways
comp_KEGG <- compareCluster(geneCluster = entrez_lists, fun = "enrichKEGG", organism = "hsa", pvalueCutoff = 0.05)
if (!is.null(comp_KEGG) && nrow(as.data.frame(comp_KEGG)) > 0) {
  
  # CHANGED: showCategory = 5 (5 per group = max 20 total on the Y-axis)
  p_comp_KEGG <- dotplot(comp_KEGG, showCategory = 5, label_format = 45, title = "Pathway Intersection (KEGG)") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  pdf(file = file.path(out_dir, "Pathway_Intersection_KEGG.pdf"), width = 10, height = 8)
  print(p_comp_KEGG)
  invisible(dev.off())
  print(p_comp_KEGG)
}