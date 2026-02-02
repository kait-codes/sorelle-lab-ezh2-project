# INSTALL ----------------------------------------------------------------------
using<-function(...) {
  
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,
                     require,
                     character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,
           require,
           character.only=TRUE)
  }
}

bioc_using<-function(...) {
  
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,
                     require,
                     character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    BiocManager::install(need)
    lapply(need,
           require,
           character.only=TRUE)
  }
}

# EDGE R -----------------------------------------------------------------------
edger_fit_genes <- function(cell_line){
  
  counts <- read.csv(paste0(data_dir, 
                            '/exp_counts/counts_', 
                            cell_line, 
                            '.csv'))
  groups <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
  dge <- DGEList(counts, group = groups) # creating DGE list
  keep <- filterByExpr(dge) #filtering out low counts
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  norm_counts <- normLibSizes(dge)
  
  png(paste0(out_dir, 
             "/mds_plots/", 
             cell_line, 
             ".png"), 
      width = 800, 
      height = 600, 
      units = "px") 
  plotMDS(norm_counts)
  dev.off()
  
  metadata_small <- read.csv(file.path(data_dir, 'samples_treatment.csv'), 
                             stringsAsFactors = TRUE)
  design <- model.matrix(~ treatment, data = metadata_small)
  rownames(design) <- colnames(norm_counts)
  
  y <- estimateDisp(norm_counts, 
                    design, 
                    robust=TRUE)
  
  fit <- glmQLFit(y, 
                  design, 
                  robust=TRUE)
  
  return(fit)
}

edger_degs <- function(fit, 
                       cell_line, 
                       comparison){
  
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
    stop('ERROR: Invalid comparison')
  }
  
  df <- as.data.frame(lrt)
  
  # converting ENSEMBL IDs to gene symbol
  genes <- df$gene_id
  annots <- select(org.Hs.eg.db, keys=genes, 
                   columns="SYMBOL", keytype="ENSEMBL")
  
  df <- merge(df, annots, 
              by.x="gene_id", by.y="ENSEMBL")
  df$SYMBOL <- ifelse(is.na(df$SYMBOL), 
                      df$gene_id, 
                      df$SYMBOL)
  
  #removing genes with PValue = 0
  df <- filter(df, PValue != 0)
  
  # annotating DEGs as no/up/down regulated
  df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > logfc_limit & df$PValue < pvalue_limit] <- "UP"
  df$diffexpressed[df$logFC < -logfc_limit & df$PValue < pvalue_limit] <- "DOWN"
  
  # saving DEG lists
  ## all DEGs
  sorted_df <- df %>% arrange(PValue)
  csv_file_path <- paste0(pipe_dir, 
                          '/degs/', 
                          cell_line, 
                          "_", 
                          comparison, 
                          ".csv")
  write.csv(sorted_df, csv_file_path, row.names = FALSE)
  
  ## significant DEGs
  sig_degs <- filter(sorted_df, diffexpressed != 'NO')
  sig_path <- paste0(pipe_dir, 
                     '/degs/', 
                     cell_line, 
                     "_", 
                     comparison, 
                     "_sig.csv")
  write.csv(sig_degs, sig_path, row.names = FALSE)
  
  ## significant up regulated DEGs
  up_degs <- filter(sorted_df, diffexpressed == 'UP')
  up_path <- paste0(pipe_dir, 
                    '/degs/', 
                    cell_line, 
                    "_", 
                    comparison, 
                    "_sig_up.csv")
  write.csv(up_degs, up_path, row.names = FALSE)
  
  ## significant down regulated DEGs
  down_degs <- filter(sorted_df, diffexpressed == 'DOWN')
  down_path <- paste0(pipe_dir, 
                      '/degs/', 
                      cell_line, 
                      "_", 
                      comparison, 
                      "_sig_down.csv")
  write.csv(down_degs, down_path, row.names = FALSE)
  
  print(paste0('Cell line: ', cell_line))
  print(paste0('Comparison: ', drugs))
  print(paste0('DEG files saved at: ', pipe_dir, '/degs'))
  
  return(df)
}



# PLOTS ------------------------------------------------------------------------
## VOLCANO PLOT -------------------------------------------
degs_volcano_plot <- function(group, 
                              comparison, 
                              label = 10, 
                              genes = NA, 
                              hover = FALSE){
  
  # defining comparison
  if (comparison == '2v1') {
    title1 <- paste0(toupper(group), ' cells')
    title2 <- 'MS177 vs DMSO'
    
  } else if (comparison == '3v1'){
    title1 <- paste0(toupper(group), ' cells')
    title2 <- 'TAZ vs DMSO'
    
  } else if (comparison == '3v2'){
    title1 <- paste0(toupper(group), ' cells')
    title2 <- 'TAZ vs MS177'
    
  } else if(comparison == 'cf5_vs_akata'){
    title1 <- paste0(toupper(group), ' treatment')
    title2 <- 'CF5 vs AKATA cells'
    
  } else {
    stop('ERROR: Invalid comparison')
    
  }
  
  print(paste0('Group: ', title1))
  print(paste0('Comparison: ', title2))
  
  # read in data
  df <- read.csv(paste0(pipe_dir, 
                        '/degs/', 
                        group, 
                        '_', 
                        comparison, 
                        '.csv'))
  
  # create a new column for labeling points based 
  # on label argument & optional gene list
  if (is.numeric(label)){
    df$delabel <- ifelse(df$SYMBOL %in% head(df[order(df$PValue), "SYMBOL"], label), 
                         df$SYMBOL, 
                         NA)
    label_title <- paste0('top', label)
    
  } else if (label == 'viral' | label == 'virus'){
    df$delabel <- ifelse(str_detect(df$SYMBOL, '^gene') == TRUE,
                         df$SYMBOL,
                         NA)
    label_title <- label
    
  } else if (label == 'custom'){
    df$delabel <- ifelse(df$SYMBOL %in% genes, 
                         df$SYMBOL, 
                         NA)
    label_title <- substitute(genes)
    
  } else{
    label <- 10
    label_title <- 'top10'
    print('ERROR: Invalid label arguement. Default top 10 genes used.')
  }

  plottitle <- paste0(title1, " in ", title2) # plot title
  
  xmax <- c(ceiling(max(abs(df$logFC))) + 1, 10) %>% 
    max() # determining x axis limits
  ymax <- c(round(max(-log10(df$PValue)), digits = -1) + 50, 100) %>% 
    max() # determining y axis limits
  
  rownames(df) <- make.unique(as.character(df$SYMBOL))

  # volcano plot
  v_plot <- ggplot(data = df, aes(x = logFC, 
                                  y = -log10(PValue), 
                                  col = diffexpressed, 
                                  label = delabel)) +
    geom_vline(xintercept = c(-logfc_limit, logfc_limit), 
               col = "gray", 
               linetype = 'dashed') +
    geom_hline(yintercept = -log10(pvalue_limit), 
               col = "gray", 
               linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    coord_cartesian(ylim = c(0, ymax), 
                    xlim = c(-xmax, xmax)) +
    labs(color = 'Severe', 
         x = expression("log"[2]*"FC"), 
         y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = seq(-xmax, xmax, 2)) +
    ggtitle(plottitle) +
    geom_label_repel(max.overlaps = Inf)

  plot_file_path <- paste0(out_dir, 
                           '/volcano_plots/', 
                           group, 
                           '_', 
                           comparison,
                           '_',
                           label_title,
                           '.png')
  
  ggsave(plot_file_path)
  print(paste0('Volcano plot saved at: ', plot_file_path))
  
  if (hover == TRUE){
    interactive_plot <- HoverLocator(v_plot)
    interactive_file_path <- paste0(out_dir, 
                                    '/volcano_plots/0_hoverlocator/', 
                                    group, 
                                    '_', 
                                    comparison, 
                                    '.html')
    
    htmlwidgets::saveWidget(interactive_plot, file = interactive_file_path)
    print(paste0('Hover volcano plot saved at: ', interactive_file_path))
  }
  
}



## VENN DIAGRAM -------------------------------------------
degs_venn_diagram <- function(comparison, 
                              group, 
                              direction){
  # defining group
  if (group == '2v1') {
    spaced_group <- 'MS177 vs DMSO'
    
  } else if (group == '3v1'){
    spaced_group <- 'TAZ vs DMSO'
    
  } else if (group == '3v2'){
    spaced_group <- 'TAZ vs MS177'
    
  } else if (group == 'cf5_vs_akata'){
    spaced_group <- 'CF5 vs AKATA'
  }
  
  # naming
  comparison_file_name <- paste(comparison, collapse = "_")
  group_sep <- unlist(strsplit(spaced_group, " ") )
  group_folder <- paste(group_sep, collapse = "")
  
  print(paste0('Comparison: ', toupper(comparison_file_name)))
  print(paste0('Group: ', spaced_group))
  
  # loading in DEG lists
  data_list <- list()
  for (i in 1:length(comparison)){
    df <- read.csv(paste0(pipe_dir, 
                          '/degs/', 
                          comparison[i], 
                          '_', 
                          group, 
                          '_sig_', 
                          direction, 
                          '.csv'))
    genes <- df$gene_id
    data_list[[length(data_list) + 1]] <- genes
  }
  
  # create & save list of shared DEGs with gene symbols included
  shared_genes <- as.data.frame(Reduce(intersect, data_list))
  colnames(shared_genes) <- 'gene_id'
  genes <- shared_genes$gene_id
  annots <- select(org.Hs.eg.db, keys=genes, 
                   columns="SYMBOL", keytype="ENSEMBL")
  shared_genes <- merge(shared_genes, annots, 
                        by.x="gene_id", by.y="ENSEMBL")
  shared_genes$SYMBOL <- ifelse(is.na(shared_genes$SYMBOL), 
                                shared_genes$gene_id, 
                                shared_genes$SYMBOL)
  
  csv_file_path <- paste0(pipe_dir, 
                          '/degs/shared/',
                          comparison_file_name, 
                          "_", 
                          direction, 
                          ".csv")
  
  write.csv(shared_genes, csv_file_path, row.names = FALSE)
  print(paste0('Shared DEGs file saved at: ', csv_file_path))
  
  # venn diagram
  venn_title <- paste0(toupper(direction), 
                       "regulated genes in ", 
                       spaced_group)
  
  ggVennDiagram(data_list, category.names = toupper(comparison)) +
    labs(title = venn_title) + 
    scale_fill_gradient(low = "blue", 
                        high = "red") + 
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(plot.background = element_rect(fill = "white"))
  
  venn_file <- paste0(out_dir, 
                      '/venn_diagrams/', 
                      group_folder, 
                      '/',
                      comparison_file_name, 
                      '_', 
                      direction, 
                      '.png')
  ggsave(venn_file)
  print(paste0('Venn diagram saved at: ', venn_file))
}


## MSigDB PATHWAY ENRICHMENT ANALYSIS DOT PLOT --------------
msigdb_pathway_enrichment_analysis <- function(group, 
                                               comparison, 
                                               collection = 'H',
                                               subcollection = NULL,
                                               padj_cutoff = 0.05,
                                               genecount_cutoff = 5){
  
  # defining comparison
  if (comparison == '2v1') {
    title1 <- paste0(toupper(group), ' cells')
    title2 <- 'MS177 vs DMSO'
    
  } else if (comparison == '3v1'){
    title1 <- paste0(toupper(group), ' cells')
    title2 <- 'TAZ vs DMSO'
    
  } else if (comparison == '3v2'){
    title1 <- paste0(toupper(group), ' cells')
    title2 <- 'TAZ vs MS177'
    
  } else if (comparison == 'cf5_vs_akata'){
    title1 <- paste0(toupper(group), ' treatment')
    title2 <- 'CF5 vs AKATA cells'
    
  } else {
    stop('ERROR: Invalid comparison')
    
  }
  
  
  # getting MSigDB collection title
  collections <- msigdbr_collections()
  if (!is.null(subcollection)){
    collection_name <- filter(collections, gs_collection == collection &
                              gs_subcollection == subcollection) %>%
      select(gs_collection_name) %>% 
      unlist() %>%
      paste(collapse = ".")
    pathway_title <- collection
  } else {
    collection_name <- filter(collections, gs_collection == collection) %>%
      select(gs_collection_name) %>% 
      unlist() %>%
      paste(collapse = ".")
    pathway_title <- paste0(collection, '_', subcollection)
  }
  if (collection_name == ""){
    stop('ERROR: Invalid MSigDB gene set collection and/or subcollection')
  }
  
  print(paste0('Group: ', title1))
  print(paste0('Comparison: ', title2))
  print(paste0('MSigDB collection: ', collection_name))
  
  # read in data & filter for human genes only
  df <- read.csv(paste0(pipe_dir, 
                        '/degs/', 
                        group, 
                        '_', 
                        comparison, 
                        '_sig.csv'))
  df_human <- filter(df, str_detect(df$gene_id, '^ENSG') == TRUE)
  deg_results_list <- split(df_human, df_human$diffexpressed)
  
  # fetch gMSigDB collection
  gene_sets_df <- msigdbr(species = 'Homo sapiens', 
                          collection = collection, 
                          subcollection = subcollection)
  gene_sets <- gene_sets_df %>% select(gs_name, ensembl_gene)
  
  # get background genes
  background_genes <- read.csv(file.path(data_dir, 'exp_counts/counts_all.csv')) %>% 
    select(1)
  
  # run clusterProfiler on each sub-dataframe
  res <- lapply(names(deg_results_list),
                function(x) enricher(gene = deg_results_list[[x]]$gene_id,
                                     TERM2GENE = gene_sets,
                                     universe = background_genes,
                                     pvalueCutoff = padj_cutoff,
                                     minGSSize = genecount_cutoff,
                                     qvalueCutoff = 0.2))
  names(res) <- names(deg_results_list)
  
  # convert the enrichResults to a dataframe with the pathways
  res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
  names(res_df) <- names(res)
  res_df <- do.call(rbind, res_df) %>% 
    mutate(minuslog10padj = -log10(p.adjust))
  # filter for significance
  res_df_sig <- res_df %>% 
    filter(p.adjust < padj_cutoff & Count > genecount_cutoff) 
  # save results
  enrich_file_path <- paste0(pipe_dir, 
                             '/pea/', 
                             group, 
                             '_', 
                             comparison,
                             '_',
                             collection,
                             '.csv')
  write.csv(res_df_sig, enrich_file_path)
  
  # visualizing 
  results_up <- res$UP
  title_up <- paste0("MSigDB ", 
                     collection_name, 
                     "\nUpregulated Genes\n", 
                     title1, 
                     " in ", 
                     title2)
  dot_up <- dotplot(results_up, 
                    showCategory = 15,
                    title = title_up)
  file_path_up <- paste0(out_dir, 
                         '/pea/', 
                         group, 
                         '_', 
                         comparison,
                         '_',
                         pathway_title,
                         'up.png')
  ggsave(file_path_up)
  
  results_down <- res$DOWN
  title_down <- paste0("MSigDB ", 
                       collection_name, 
                       "\nDownregulated Genes\n", 
                       title1, 
                       " in ", 
                       title2)
  dot_down <- dotplot(results_down, 
                      showCategory = 15,
                      title = title_down)
  file_path_down <- paste0(out_dir, 
                           '/pea/', 
                           group, 
                           '_', 
                           comparison,
                           '_',
                           pathway_title,
                           'down.png')
  ggsave(file_path_down)
  
  
  print(paste0('PEA plots saved at: ', out_dir, '/pea/'))
  
  ## OPTIONAL: uncomment to allow for further plot making
  #return(res)
}







