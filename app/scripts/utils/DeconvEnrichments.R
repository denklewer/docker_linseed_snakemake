library(Seurat)
library(ggplot2)

plotDeconvEnrichment <- function(SC_PATH,SC_DATASET,SAVE_PATH,PATH,DT,annotation_column = "seurat_clusters",
                                 assay = "RNA", top_genes = 100, calculate_limma = T,
                                 filter_ = T, filter_limit = 0.01, render_=F) {
  dir.create(file.path(SAVE_PATH,DT),showWarnings = F, recursive = T)
  
    folders_ <- basename(list.dirs(PATH,recursive = F))
  folders_ <- folders_[grepl("^ct",folders_)]
  cts_ <- as.numeric(gsub("^ct","",folders_))
  min_ct <- min(cts_)
  max_ct <- max(cts_)
  if (!render_) {
  sc_data <- readRDS(SC_PATH)
  expData <- get(assay,sc_data@assays)@data
  reductions <- sc_data@reductions
  metadata <- sc_data@meta.data
  rm(sc_data)
  gc()
  
  for (cell_type in min_ct:max_ct) {
    
    LOAD.PATH <- file.path(PATH,paste0("ct",cell_type),DT,"best")
    SAVE.PATH <- file.path(SAVE_PATH,DT,paste0("ct",cell_type))
    dir.create(SAVE.PATH,showWarnings = F, recursive = T)
    
    basis_ <- data.frame(fread(file.path(LOAD.PATH,"basis_column.tsv"))) %>%
      dplyr::rename(gene = "V1")
    avg_props <- t(basis_[1,-(1:(cell_type+1))])
    basis_ <- basis_[-1,]
    
    if (filter_) {
      avg_props <- avg_props[avg_props>filter_limit,]
      cell_types <- length(avg_props)
      basis_ <- basis_[,c("gene",paste0("FC_",names(avg_props)),
                          names(avg_props))]
    }
    cols_ <- colnames(basis_)[grepl("^Cell.*",colnames(basis_))]
    if (calculate_limma) {
      limma_basis <- matrix(0,nrow=nrow(basis_),ncol=cell_types)
      colnames(limma_basis) <- paste0("Limma_FC_",cols_)
      for (col_ in cols_) {
        expr <- basis_[,cols_]
        condition <- as.integer(col_ != colnames(expr))
        design <- model.matrix(~ condition)
        fit <- lmFit(expr, design)
        fit <- eBayes(fit)
        stats <- topTable(fit, number = nrow(basis_))
        limma_basis[,paste0("Limma_FC_",col_)] <- -stats[rownames(basis_),"logFC"]
      }
      basis_ <- cbind(limma_basis,basis_)  
    }
    fc_column <- "FC_"
    if (calculate_limma) {
      fc_column <- "Limma_FC_"
    }
    
    for (col_ in cols_) {
      genes <- (basis_ %>% arrange(desc(.data[[paste0(fc_column,col_)]])))[1:top_genes,"gene"]
      genes <- genes[tolower(genes) %in% tolower(rownames(expData))]
      write.table(data.frame(genes), file.path(SAVE.PATH,paste0("markers_",tolower(col_),".txt")),
                  row.names = F,col.names = F,quote = F)
      k <- length(genes)
      av <- numeric(ncol(expData))
      
      zz <- which(tolower(rownames(expData)) %in% tolower(genes))
      geneExp <- as.matrix(expData[zz, ])
      avgExpData <- apply(geneExp,2,mean)
      
      if ("tsne" %in% names(reductions)) {
        
        dir.create(file.path(SAVE.PATH,"SCN",SC_DATASET,"tSNE"),showWarnings = F, recursive = T)
        
        toPlot <- data.frame(tSNE_1 = reductions$tsne@cell.embeddings[,1], 
                             tSNE_2 = reductions$tsne@cell.embeddings[,2],
                             avg_exp = avgExpData,
                             cluster = metadata[,annotation_column])
        
        centers <- toPlot %>% 
          dplyr::group_by_at("cluster") %>% 
          dplyr::summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) 
        
        pl <- ggplot(data=toPlot, aes_string(x="tSNE_1", y="tSNE_2",color="avg_exp")) + geom_point() +
          theme_classic() + theme(legend.text=element_text(size=5)) +
          expand_limits(x=0, y=0) +
          ggtitle(paste0(col_," - ",sprintf('%.2f%%',avg_props[col_]*100))) + coord_fixed(2/3) +
          scale_color_gradientn(colors=c("chartreuse","darkorchid3")) +
          geom_text(data = centers, aes_string(label="cluster"), size = 3, color = "black")
        
        png(file.path(SAVE.PATH,"SCN",SC_DATASET,"tSNE",paste0("SC_",tolower(col_),".png")))
        print(pl)
        dev.off()
      }
      
      if ("umap" %in% names(reductions)) {
        dir.create(file.path(SAVE.PATH,"SCN",SC_DATASET,"UMAP"),showWarnings = F, recursive = T)
        toPlot <- data.frame(UMAP_1 = reductions$umap@cell.embeddings[,1], 
                             UMAP_2 = reductions$umap@cell.embeddings[,2],
                             avg_exp = avgExpData,
                             cluster = metadata[,annotation_column])
        
        centers <- toPlot %>% 
          dplyr::group_by_at("cluster") %>% 
          dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) 
        
        pl <- ggplot(data=toPlot, aes_string(x="UMAP_1", y="UMAP_2",color="avg_exp")) + geom_point() +
          theme_classic() + theme(legend.text=element_text(size=5)) +
          expand_limits(x=0, y=0) +
          ggtitle(paste0(col_," - ",sprintf('%.2f%%',avg_props[col_]*100))) + coord_fixed(2/3) +
          scale_color_gradientn(colors=c("chartreuse","darkorchid3")) +
          geom_text(data = centers, aes_string(label="cluster"), size = 3, color = "black")
        
        png(file.path(SAVE.PATH,"SCN",SC_DATASET,"UMAP",paste0("SC_",tolower(col_),".png")))
        print(pl)
        dev.off()
      }
    }
  } 
  }
  #"/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/Docker/docker_linseed_snakemake/app/scripts/utils/reportSCReference.Rmd",
  suppressMessages(rmarkdown::render(#"/app/scripts/utils/reportSCReference.Rmd",
    "/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/Docker/docker_linseed_snakemake/app/scripts/utils/reportSCReference.Rmd",
                    output_dir = file.path(SAVE_PATH,DT), output_file = paste0(DT,".html"),
                    params = list(sc_dataset = SC_DATASET,
                                  res_data = SC_PATH,
                                  res_type = "SCN",
                                  res_path = file.path(SAVE_PATH,DT),
                                  min_ct = min_ct,
                                  max_ct = max_ct),quiet = T))
}

# Example
# SC_PATH <- "/Volumes/martyomov/Active/IndividualBackUps/alexeys/deconvolution/single_cell/GSE159115_seurat.rds"
# SC_DATASET <- "GSE159115"
# SAVE_PATH <- "/Volumes/martyomov/Active/IndividualBackUps/alexeys/deconvolution/data/10_2022/reports/KIRC_RPLS_CODING_MAD1_KNN_MOVEMENT_LIMITS/GSE159115"
# PATH <- "/Volumes/martyomov/Active/IndividualBackUps/alexeys/deconvolution/data/10_2022/results/KIRC_RPLS_CODING_MAD1_KNN_MOVEMENT_LIMITS"
# DT <- "20221007"
# plotDeconvEnrichment(SC_PATH,SC_DATASET,SAVE_PATH,PATH,DT,assay = "SCT",render=T)
