# Functions for analyzing RNAseq of EBV+ samples
# Author: Kaitlyn Tremble
# Last updated: 2026-02-23

# PACKAGE INSTALL & LOADING ----------------------------------------------------
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



nice_title <- function(group, comparison){
  # defining comparison for nice title
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
  
  title_list <- list(title1, title2)
  return(title_list)
}

# EDGE R -----------------------------------------------------------------------
# see BioConductor edgeR user manual 

# first fit gene counts to a model
edger_fit_genes <- function(cell_line){
  
  # read in expected counts file
  counts <- read.csv(paste0(data_dir, 
                            '/exp_counts/counts_', 
                            cell_line, 
                            '.csv'))
  
  # see manual 2.6
  ## change if number of samples in each group is different
  groups <- rep(1:3, each = 4)
  dge <- DGEList(counts, group = groups) #creating DGE list
  
  # see manual 2.7
  keep <- filterByExpr(dge) #filtering out low counts
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  # normalizing by library size
  # see manual 2.8.3
  norm_counts <- normLibSizes(dge)
  
  # MDS plot of samples
  png(paste0(out_dir, 
             "/samples/", 
             cell_line, 
             "_mds.png"), 
      width = 800, 
      height = 600, 
      units = "px") 
  my_colors <- c('blue', 'red', 'purple')
  color_vec <- rep(my_colors, each = 4)
  plotMDS(norm_counts, col= color_vec)
  legend("topright", col=my_colors, legend=c("DMSO", "MS177", "TAZ"),
         text.col=my_colors)
  dev.off()
  
  # creating design matrix
  # note: only include factors you will use (in this case treatment)
  # my matrix doesn't include sample labels since I am using it for multiple
  # cell lines/sample labels
  ## I referenced this tutorial for this:
  ## https://gtpb.github.io/ADER18S/pages/tutorial_complex
  metadata_small <- read.csv(file.path(data_dir, 'samples_treatment.csv'),
                             stringsAsFactors = TRUE)
  design <- model.matrix(~ treatment, data = metadata_small)
  rownames(design) <- colnames(norm_counts) #ensure proper grouping here

  # using the negative binomial GLM framework to estimate gene dispersions
  # see manual 2.11.2
  y <- estimateDisp(norm_counts,
                    design,
                    robust=TRUE)

  fit <- glmQLFit(y,
                  design,
                  robust=TRUE)

  return(fit)
}

# use output of above function to identify DEGs
edger_degs <- function(fit, 
                       cell_line, 
                       comparison){
  
  # defining drug comparison (based on design matrix order so be careful)
  # see manual 2.11.3
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
  
  # removing genes with PValue = 0
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
  print('---------------------------------------------------------------------')
  
  # OPTIONAL: use DEGs directly from this function for further analysis
  # I code in the intermediates as saved .csv files so I don't do this
  #return(df)
}



# PLOTS ------------------------------------------------------------------------
## VOLCANO PLOT -------------------------------------------
# I referenced the BioStat Squid tutorial on volcano plots

## can label viral genes using label = 'viral'
## or label custom gene list using label = 'custom' and gene = your_list
## can also use hover = TRUE to make interactive plot

degs_volcano_plot <- function(group, 
                              comparison, 
                              label = 10, 
                              genes = NA, 
                              hover = FALSE){
  
  title_list <- nice_title(group, comparison)
  
  print(paste0('Group: ', title_list[1]))
  print(paste0('Comparison: ', title_list[2]))
  
  # read in data (make sure your naming convention is consistent)
  df <- read.csv(paste0(pipe_dir, 
                        '/degs/', 
                        group, 
                        '_', 
                        comparison, 
                        '.csv'))
  
  # create a new column for labeling points based 
  # on label argument & optional gene list
  if (is.numeric(label)){
    df$delabel <- ifelse(df$SYMBOL %in% head(df[order(df$PValue), 
                                                "SYMBOL"], label), 
                         df$SYMBOL, 
                         NA)
    label_title <- paste0('top', label)
    
  } else if (label == 'viral' | label == 'virus'){
    df$delabel <- ifelse(str_detect(df$SYMBOL, '^ENSG') == FALSE,
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

  plottitle <- paste0(title_list[1], " in ", title_list[2]) #plot title
  
  xmax <- c(ceiling(max(abs(df$logFC))) + 1, 10) %>% 
    max() #determining x axis limits
  ymax <- c(round(max(-log10(df$PValue)), digits = -1) + 50, 100) %>% 
    max() #determining y axis limits
  
  rownames(df) <- make.unique(as.character(df$SYMBOL)) #for Hover plot

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
                       labels = c("Downregulated", 
                                  "Not significant", 
                                  "Upregulated")) +
    coord_cartesian(ylim = c(0, ymax), 
                    xlim = c(-xmax, xmax)) +
    labs(color = 'Severe', 
         x = expression("log"[2]*"FC"), 
         y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = seq(-xmax, xmax, 2)) +
    ggtitle(plottitle) +
    geom_label_repel(max.overlaps = Inf) +
    theme_bw()

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
  
  # Hover volcano plot
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
  
  print('---------------------------------------------------------------------')
  ## OPTIONAL: uncomment to continue working with the plot in R
  # return(v_plot)
  
}


## VENN DIAGRAM -------------------------------------------
degs_venn_diagram <- function(comparison, 
                              group, 
                              direction = c('up', 'down')){
  
  # defining group for nice title
  if (group == '2v1') {
    spaced_group <- 'MS177 vs DMSO'
    
  } else if (group == '3v1'){
    spaced_group <- 'TAZ vs DMSO'
    
  } else if (group == '3v2'){
    spaced_group <- 'TAZ vs MS177'
    
  } else if (group == 'cf5_vs_akata'){
    spaced_group <- 'CF5 vs AKATA'
  }
  
  # naming file paths
  comparison_file_name <- paste(comparison, collapse = "_")
  group_sep <- unlist(strsplit(spaced_group, " ") )
  group_folder <- paste(group_sep, collapse = "")
  
  print(paste0('Comparison: ', toupper(comparison_file_name)))
  print(paste0('Group: ', spaced_group))
  
  # loading in DEG lists (significant only)
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
                          group,
                          "_",
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
  print('---------------------------------------------------------------------')
}


## ANALYSIS USING MSigDB -------------------------------------------------------
# can change collection & subcollection
# use msigdbr_collections() to view options

title_msigdb <- function(collection = 'H',
                         subcollection = NULL){
  
  # getting MSigDB collection title
  collections <- msigdbr_collections()
  if (!is.null(subcollection)){
    collection_name <- filter(collections, gs_collection == collection &
                                gs_subcollection == subcollection) %>%
      dplyr::select(gs_collection_name) %>% 
      unlist() %>%
      paste(collapse = ".")
    pathway_title <- paste0(collection, '_', subcollection)
  } else {
    collection_name <- filter(collections, gs_collection == collection) %>%
      dplyr::select(gs_collection_name) %>% 
      unlist() %>%
      paste(collapse = ".")
    pathway_title <- collection
  }
  if (collection_name == ""){
    stop('ERROR: Invalid MSigDB gene set collection and/or subcollection')
  }
  
  msigdb_list <- list(collection_name, pathway_title)
  return(msigdb_list)
}


fetch_msigdb <- function(collection = 'H',
                         subcollection = NULL){
  
  gene_sets_df <- msigdbr(species = 'Homo sapiens', 
                          collection = collection, 
                          subcollection = subcollection)
  
  return(gene_sets_df)
}


## PATHWAY ENRICHMENT ANALYSIS DOT PLOT ---------------------
## I referenced the BioStat Squid tutorial:
## Pathway Enrichment Analysis with clusterProfiler
msigdb_pathway_enrichment_analysis <- function(group, 
                                               comparison, 
                                               collection = 'H',
                                               subcollection = NULL,
                                               padj_cutoff = 0.05,
                                               genecount_cutoff = 5){
  
  # nice titles
  title_list <- nice_title(group, comparison)
  
  # getting MSigDB collection title
  msigdb_list <- title_msigdb(collection, subcollection)
  
  # fetch MSigDB collection
  gene_sets_df <- fetch_msigdb(collection, subcollection)
  gene_sets <- gene_sets_df %>% dplyr::select(gs_name, ensembl_gene)
  
  print(paste0('Group: ', title_list[1]))
  print(paste0('Comparison: ', title_list[2]))
  print(paste0('MSigDB collection: ', msigdb_list[1]))
  
  # read in DEG list & split by up/down regulated
  df <- read.csv(paste0(pipe_dir, 
                        '/degs/', 
                        group, 
                        '_', 
                        comparison, 
                        '_sig.csv'))
  deg_results_list <- split(df, df$diffexpressed)
  
  # get background genes from original exp counts file
  background_genes <- read.csv(file.path(data_dir, 'exp_counts/counts_all.csv')) %>% 
    dplyr::select(1)
  
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
                             msigdb_list[2],
                             '.csv')
  write.csv(res_df_sig, enrich_file_path)
  print(paste0('PEA results saved at: ', enrich_file_path))
  
  # visualizing results for up & down regulated genes
  results_up <- res$UP
  title_up <- paste0("MSigDB ", 
                      msigdb_list[1], 
                      "\nUpregulated Genes\n", 
                      title_list[1], 
                      " in ", 
                      title_list[2])
  dot_up <- dotplot(results_up, 
                    showCategory = 15,
                    title = title_up)
  file_path_up <- paste0(out_dir, 
                          '/pea/', 
                          group, 
                          '_', 
                          comparison,
                          '_',
                         msigdb_list[2],
                          '_up.png')
  
  if (length(dot_up[["data"]][["ID"]]) == 0) {
    print("No UPregulated pathways")
    file.remove(file_path_up)
  } else {
    ggsave(file_path_up)
  }

  
  results_down <- res$DOWN
  title_down <- paste0("MSigDB ", 
                       msigdb_list[1], 
                        "\nDownregulated Genes\n", 
                        title_list[1], 
                        " in ", 
                        title_list[2])
  dot_down <- dotplot(results_down, 
                      showCategory = 15,
                      title = title_down)
  file_path_down <- paste0(out_dir, 
                            '/pea/', 
                            group, 
                            '_', 
                            comparison,
                            '_',
                           msigdb_list[2],
                            '_down.png')
    
  if (length(dot_down[["data"]][["ID"]]) == 0) {
    print("No DOWNregulated pathways")
    file.remove(file_path_down)
  } else {
    ggsave(file_path_down)
  }
  
  
  print(paste0('PEA plots saved at: ', out_dir, '/pea/'))
  print('---------------------------------------------------------------------')
  
  ## OPTIONAL: uncomment to allow for further plot making besides dot plot
  #return(res)

  ## or uncomment to further edit the dot plots in R
  #plot_list <- list(dot_up, dot_down)
  #return(plot_list)
}



## GENE SET ENRICHMENT ANALYSIS -----------------------------


gene_ranking <- function(group, 
                         comparison){
  
  # read in DEG list
  df <- read.csv(paste0(pipe_dir, 
                        '/degs/', 
                        group, 
                        '_', 
                        comparison, 
                        '.csv'))
  
  # filter out viral genes & fix gene ids
  df <- df %>% dplyr::filter(str_detect(gene_id, "^ENSG"))
  df$gene_id <- make.unique(as.character(df$gene_id))
  
  # rank genes by logFC & p-value
  rankings <- sign(df$logFC)*(-log10(df$PValue))
  names(rankings) <- df$gene_id
  
  # fix infinite values caused by small p-values
  max_ranking <- max(rankings[is.finite(rankings)])
  min_ranking <- min(rankings[is.finite(rankings)])
  rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
  rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
  
  rankings <- sort(rankings, decreasing = TRUE) #sort genes by ranking
  # plot(rankings) #uncomment to check rank sort
  
  ## uncomment & adjust to check rankings of specific genes
  ## (currently first 50 by alphabetical order)
  # ranked_genes_plot <- ggplot(data.frame(gene_symbol = names(rankings)[1:50], 
  #                                   ranks = rankings[1:50]), 
  #                        aes(gene_symbol, ranks)) + 
  #   geom_point() +
  #   theme_classic() + 
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # ranked_genes_plot
  
  return(rankings)
}


msigdb_gsea <- function(group, 
                        comparison, 
                        collection = 'H',
                        subcollection = NULL){
  
  # nice title
  title_list <- nice_title(group, comparison)
  
  # getting MSigDB collection title
  msigdb_list <- title_msigdb(collection, subcollection)
  
  # fetch MSigDB collection
  gene_sets_df <- fetch_msigdb(collection, subcollection)
  gene_sets <- gene_sets_df %>% dplyr::select(gs_name, ensembl_gene)
  
  print(paste0('Group: ', title_list[1]))
  print(paste0('Comparison: ', title_list[2]))
  print(paste0('MSigDB collection: ', msigdb_list[1]))
  
  # rank genes
  rankings <- gene_ranking(group, comparison)
  
  # run fgsea
  GSEAres <- GSEA(rankings, 
                  TERM2GENE = gene_sets)
  
  # save GSEA results
  gsea_results_file <- paste0(pipe_dir, 
                              '/gsea/', 
                              group, 
                              '_', 
                              comparison,
                              '_',
                              msigdb_list[2],
                              '.RDS')
  
  saveRDS(GSEAres, gsea_results_file)
  print(paste0('GSEA results saved at: ', gsea_results_file))
  print('---------------------------------------------------------------------')
  
  return(GSEAres)
}


## for this I referenced the BioStat Squid tutorial:
## Easy Gene Set Enrichment Analysis in R with fgsea()
msigdb_gsea_main_plot <- function(group, 
                           comparison, 
                           collection = 'H',
                           subcollection = NULL){
  
  # nice title
  title_list <- nice_title(group, comparison)
  
  # getting MSigDB collection title
  msigdb_list <- title_msigdb(collection, subcollection)
  
  # fetch MSigDB collection
  gene_sets_df <- fetch_msigdb(collection, subcollection)
  gene_sets <- gene_sets_df %>%
    split(x = .$ensembl_gene, f = .$gs_name)
  
  # rank genes
  rankings <- gene_ranking(group, comparison)
  
  # run fgsea
  GSEAdf <- fgsea(pathways = gene_sets, # list of gene sets to check
                   stats = rankings,
                   scoreType = 'std',
                   minSize = 10,
                   maxSize = 500,
                   nproc = 1) # for parallelisation
  
  # select only independent pathways, removing redundancies/similar pathways
  collapsedPathways <- collapsePathways(GSEAdf[order(pval)][pval < 0.05], 
                                        gene_sets, rankings)
  mainPathways <- GSEAdf[pathway %in% 
                            collapsedPathways$mainPathways][order(-NES), pathway]
  
  gsea_plot_file <- paste0(out_dir, 
                           '/gsea/', 
                           group, 
                           '_', 
                           comparison,
                           '_',
                           msigdb_list[2],
                           '_main.pdf')
  
  # plot top pathways
  p <- plotGseaTable(gene_sets[head(mainPathways, n=20)], 
                     rankings, GSEAdf, gseaParam = 0.5) +
  pdf(file = gsea_plot_file, width = 20, height = 12)
  print(p)
  dev.off()
  
  print(paste0('GSEA plot saved at: ', gsea_plot_file))
  print('---------------------------------------------------------------------')

}


msigdb_gsea_enrich_plot <- function(name, 
                                    gsea_result,
                                    pathway_num){
  
  filename <- paste0(out_dir, 
                     '/gsea/',
                     name,
                     '_',
                     gsea_result$Description[pathway_num],
                     '.png')
  
  p <- gseaplot2(gsea_result,
                 geneSetID = pathway_num,
                 title = gsea_result$Description[pathway_num],
                 base_size = 20)
  
  png(file = filename, width = 800, height = 800)
  print(p)
  dev.off()

  print(paste0('GSEA enrich plot saved at: ', filename))
  print('---------------------------------------------------------------------')
  
}

