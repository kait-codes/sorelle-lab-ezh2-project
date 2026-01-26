# Initial package install
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


edger_fit_genes <- function(cell_line){
  counts <- read.csv(paste0(pipe_dir, '/exp_counts/counts_', cell_line, '.csv'))
  dge <- DGEList(counts) # creating DGE list
  keep <- filterByExpr(dge) #filtering out low counts
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  norm_counts <- normLibSizes(dge)
  
  png(paste0(out_dir, "/mds_plots/", cell_line, "_mds.png"), 
      width = 800, height = 600, units = "px") 
  plotMDS(norm_counts)
  dev.off()
  
  metadata_small <- read.csv(file.path(data_dir, 'samples2.csv'), stringsAsFactors = TRUE)
  design <- model.matrix(~ treatment, data = metadata_small)
  rownames(design) <- colnames(norm_counts)
  y <- estimateDisp(norm_counts, design, robust=TRUE)
  plotBCV(y)
  fit <- glmQLFit(y, design, robust=TRUE)
  
  return(fit)
}

edger_volcano_plot <- function(fit, cell_line, comparison){
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
  
  # annotating DEGs as no/up/down regulated
  df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > 0.6 & df$PValue < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -0.6 & df$PValue < 0.05] <- "DOWN"
  
  # create a new column for the names of the top 30 DEGs
  df$delabel <- ifelse(df$SYMBOL %in% head(df[order(df$PValue), "SYMBOL"], 30), df$SYMBOL, NA)
  
  sorted_df <- df %>% arrange(PValue)
  csv_file_path <- paste0(pipe_dir, '/degs/', cell_line, "_", comparison, "_deg.csv")
  write.csv(sorted_df, csv_file_path)
  

  plottitle <- paste0(toupper(cell_line), " cells in ", drugs)
  
  # volcano plot
  ggplot(data = df, aes(x = logFC, y = -log10(PValue), col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 2) + 
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # set the colors  
                       labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels 
    labs(color = 'Severe', # legend title 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # customize breaks in the x axis
    ggtitle(plottitle) + # plot title
    geom_text() ##FIXME idk how to make this look better
  
  plot_file_path <- paste0(out_dir, '/volcano_plots/', cell_line, '_', comparison, '_volcano.png')
  ggsave(plot_file_path)
  
  print(paste0('Comparison: ', drugs))
  print(paste0('DEG file saved at: ', csv_file_path))
  print(paste0('Volcano plot saved at: ', plot_file_path))
}