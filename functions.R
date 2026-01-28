# INSTALL ----------------------------------------------------------------------
using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

bioc_using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    BiocManager::install(need)
    lapply(need,require,character.only=TRUE)
  }
}

# EDGE R -----------------------------------------------------------------------
edger_fit_genes <- function(cell_line){
  counts <- read.csv(paste0(data_dir, '/exp_counts/counts_', cell_line, '.csv'))
  groups <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
  dge <- DGEList(counts, group = groups) # creating DGE list
  keep <- filterByExpr(dge) #filtering out low counts
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  norm_counts <- normLibSizes(dge)
  
  png(paste0(out_dir, "/mds_plots/", cell_line, ".png"), 
      width = 800, height = 600, units = "px") 
  plotMDS(norm_counts)
  dev.off()
  
  metadata_small <- read.csv(file.path(data_dir, 'samples_treatment.csv'), stringsAsFactors = TRUE)
  design <- model.matrix(~ treatment, data = metadata_small)
  rownames(design) <- colnames(norm_counts)
  y <- estimateDisp(norm_counts, design, robust=TRUE)
  plotBCV(y)
  fit <- glmQLFit(y, design, robust=TRUE)
  
  return(fit)
}

edger_degs <- function(fit, cell_line, comparison){
  # defining drug comparison
  if (comparison == '2v1') {
    lrt <- glmLRT(fit, coef=2)
    drugs <- 'MS177 vs DMSO'
  } else if (comparison == '3v1'){
    lrt <- glmLRT(fit, coef=3)
    drugs <- 'TAZ vs DMSO'
  } else if (comparison == '3v2'){
    lrt <- glmLRT(fit, contrast=c(0,-1,1))
    drugs <- 'TAZ vs MS177'
  } else {
    print('ERROR: Invalid comparison')
  }
  
  ## should put in condition so function stops here without valid comparison
  df <- as.data.frame(lrt)
  
  # converting ENSEMBL IDs to gene symbol
  genes <- df$gene_id
  annots <- select(org.Hs.eg.db, keys=genes, 
                   columns="SYMBOL", keytype="ENSEMBL")
  
  df <- merge(df, annots, by.x="gene_id", by.y="ENSEMBL")
  df$SYMBOL <- ifelse(is.na(df$SYMBOL), df$gene_id, df$SYMBOL)
  
  #removing genes with PValue = 0
  df <- filter(df, PValue != 0)
  
  # annotating DEGs as no/up/down regulated
  df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > logfc_limit & df$PValue < pvalue_limit] <- "UP"
  df$diffexpressed[df$logFC < -logfc_limit & df$PValue < pvalue_limit] <- "DOWN"
  
  # saving DEG lists
  ## all DEGs
  sorted_df <- df %>% arrange(PValue)
  csv_file_path <- paste0(pipe_dir, '/degs/', cell_line, "_", comparison, ".csv")
  write.csv(sorted_df, csv_file_path, row.names = FALSE)
  ## significant DEGs
  sig_degs <- filter(sorted_df, diffexpressed != 'NO')
  sig_path <- paste0(pipe_dir, '/degs/', cell_line, "_", comparison, "_sig.csv")
  write.csv(sig_degs, sig_path, row.names = FALSE)
  ## significant up regulated DEGs
  up_degs <- filter(sorted_df, diffexpressed == 'UP')
  up_path <- paste0(pipe_dir, '/degs/', cell_line, "_", comparison, "_sig_up.csv")
  write.csv(up_degs, up_path, row.names = FALSE)
  ## significant down regulated DEGs
  down_degs <- filter(sorted_df, diffexpressed == 'DOWN')
  down_path <- paste0(pipe_dir, '/degs/', cell_line, "_", comparison, "_sig_down.csv")
  write.csv(down_degs, down_path, row.names = FALSE)
  
  print(paste0('Cell line: ', cell_line))
  print(paste0('Comparison: ', drugs))
  print(paste0('DEG files saved at: ', pipe_dir, '/degs'))
  
  return(df)
}

# PLOTS ------------------------------------------------------------------------
degs_volcano_plot <- function(cell_line, comparison){
  # defining drug comparison
  if (comparison == '2v1') {
    lrt <- glmLRT(fit, coef=2)
    drugs <- 'MS177 vs DMSO'
  } else if (comparison == '3v1'){
    lrt <- glmLRT(fit, coef=3)
    drugs <- 'TAZ vs DMSO'
  } else if (comparison == '3v2'){
    lrt <- glmLRT(fit, contrast=c(0,-1,1))
    drugs <- 'TAZ vs MS177'
  } else {
    print('ERROR: Invalid comparison')
  }
  
  print(paste0('Cell line: ', cell_line))
  print(paste0('Comparison: ', drugs))
  
  # read in data
  df <- read.csv(paste0(pipe_dir, '/degs/', cell_line, '_', comparison, '.csv'))
  
  # create a new column for the names of the top 10 DEGs
  df$delabel <- ifelse(df$SYMBOL %in% head(df[order(df$PValue), "SYMBOL"], 10), df$SYMBOL, NA)

  plottitle <- paste0(toupper(cell_line), " cells in ", drugs)
  options <- c(round(max(-log10(df$PValue)), digits = -1) + 50, 100)
  ymax <- max(options)

  # volcano plot
  ggplot(data = df, aes(x = logFC, y = -log10(PValue), col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = c(-logfc_limit, logfc_limit), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(pvalue_limit), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # set the colors
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    coord_cartesian(ylim = c(0, ymax), xlim = c(-10, 10)) +
    labs(color = 'Severe', # legend title
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # customize breaks in the x axis
    ggtitle(plottitle) + # plot title
    geom_label_repel(max.overlaps = Inf)

  plot_file_path <- paste0(out_dir, '/volcano_plots/', cell_line, '_', 
                           comparison, '.png')
  ggsave(plot_file_path)
  print(paste0('Volcano plot saved at: ', plot_file_path))
}


degs_venn_diagram <- function(compare_cells, treatment, direction){
  # defining treatment comparison
  if (treatment == '2v1') {
    drugs <- 'MS177 vs DMSO'
  } else if (treatment == '3v1'){
    drugs <- 'TAZ vs DMSO'
  } else if (treatment == '3v2'){
    drugs <- 'TAZ vs MS177'
  }
  
  # naming
  cells_file_name <- paste(compare_cells, collapse = "_")
  drugs_sep <- unlist(strsplit(drugs, " ") )
  drugs_folder <- paste(drugs_sep, collapse = "")
  print(paste0('Cell lines: ', toupper(cells_file_name)))
  print(paste0('Treatment comparison: ', drugs))
  
  # loading in DEG lists
  data_list <- list()
  for (i in 1:length(compare_cells)){
    df <- read.csv(paste0(pipe_dir, '/degs/', compare_cells[i], '_', treatment, 
                          '_sig_', direction, '.csv'))
    genes <- df$gene_id
    data_list[[length(data_list) + 1]] <- genes
  }
  
  # create & save list of shared DEGs
  shared_genes <- as.data.frame(Reduce(intersect, data_list))
  colnames(shared_genes) <- paste0(drugs_folder, '_', cells_file_name, '_', direction)
  csv_file_path <- paste0(pipe_dir, '/degs/shared/',
                          cells_file_name, "_", direction, ".csv")
  write.csv(shared_genes, csv_file_path, row.names = FALSE)
  print(paste0('Shared DEGs file save at: ', csv_file_path))
  
  # venn diagram
  venn_title <- paste0(toupper(direction), "regulated genes in ", drugs)
  ggVennDiagram(data_list, category.names = toupper(compare_cells)) +
    labs(title = venn_title) + 
    scale_fill_gradient(low = "blue", high = "red") + 
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(plot.background = element_rect(fill = "white"))

  venn_file <- paste0(out_dir, '/venn_diagrams/', drugs_folder, '/',
                      cells_file_name, '_', direction, '.png')
  ggsave(venn_file)
  print(paste0('Venn diagram saved at: ', venn_file))
}
