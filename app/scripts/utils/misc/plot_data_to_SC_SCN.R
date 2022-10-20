library(rhdf5)
library(rjson)
library(jsonlite)
library(Matrix)
library(dplyr)
library(data.table)
library(ggplot2)
library(optparse)
library(limma)

option_list = list(
  make_option(c("--dataset"), type="character", default=NULL),
  make_option(c("--results"), type="character", default=NULL),
  make_option(c("--save"), type="character", default=getwd()),
  make_option(c("--min_ct"), type="numeric", default=NULL),
  make_option(c("--max_ct"), type="numeric", default=NULL),
  make_option(c("--top_genes"), type="numeric", default=100),
  make_option(c("--dt"), type="character", default=NULL),
  make_option(c("--column"), type="character", default="Cluster")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


SC_PATH <- opt$dataset #file.path("/Users/aladyeva.e/Downloads/GSE131907")
SC_DATASET <- basename(SC_PATH) #"GSE131907"
RES_PATH <- opt$results #file.path("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/Snakemake/reports",RESULTS_DIR)
RESULTS_DIR <- strsplit(RES_PATH,"/")[[1]] #basename(RES_PATH)#"LUAD_F_HK_CC_MAD32_GENES_2000"
RESULTS_DIR <- RESULTS_DIR[length(RESULTS_DIR)-1]
SAVE.PATH <- opt$save

min_ct <- opt$min_ct
max_ct <- opt$max_ct
top_genes <- opt$top_genes

DT <- opt$dt
if (is.null(opt$dt)){
  DT <- list.files(file.path(RES_PATH,paste0("ct",opt$min_ct)))[1]
}
print(DT)

result <- jsonlite::fromJSON(file.path(SC_PATH,"exp_data.json"))
total_counts <- result$totalCounts
data_ <- h5read(file.path(SC_PATH,"data.h5"),"X")
dims_ <- h5readAttributes(file = file.path(SC_PATH,"data.h5"), 
                          name = "X")$shape
expData <- sparseMatrix(p = data_$indptr, j = data_$indices, x = as.numeric(data_$data), 
                        dimnames = list(result$features, result$barcodes),
                        dims=c(dims_[2],dims_[1]), repr="R", index1=FALSE)
plotData <- (jsonlite::fromJSON(file.path(SC_PATH,"plot_data.json"))$data) %>% as.data.frame

dir.create(file.path(opt$save,DT),showWarnings = F, recursive = T)

for (cell_types in min_ct:max_ct) {
  
  LOAD.PATH <- file.path(RES_PATH,paste0("ct",cell_types),DT,"best")
  SAVE.PATH <- file.path(opt$save,DT,paste0("ct",cell_types))
  
  basis_ <- data.frame(fread(file.path(LOAD.PATH,"basis_column.tsv"))) %>%
    dplyr::rename(gene = "V1")
  avg_props <- t(basis_[1,-(1:(cell_types+1))])
  basis_ <- basis_[-1,]
  
  dir.create(file.path(SAVE.PATH,"SCN",SC_DATASET,"UMAP"),showWarnings = F, recursive = T)
  dir.create(file.path(SAVE.PATH,"SCN",SC_DATASET,"tSNE"),showWarnings = F, recursive = T)
  
  for (ct in 1:cell_types) {
    expr <- basis_[,grepl("^Cell.*",colnames(basis_))]
    condition <- as.integer(paste0("Cell_type_",ct) != colnames(expr))
    design <- model.matrix(~ condition)
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    stats <- topTable(fit, number = nrow(basis_))
    basis_[,paste0("FC_Cell_type_",ct)] <- -stats[rownames(basis_),"logFC"]
  }
  
  for (ct in 1:cell_types) {
    genes <- (basis_ %>% arrange(desc(.data[[paste0("FC_Cell_type_",ct)]])))[1:top_genes,"gene"]
    genes <- genes[tolower(genes) %in% tolower(rownames(expData))]
    write.table( data.frame(genes), file.path(SAVE.PATH,paste0("markers_cell_type_",ct,".txt")),
                 row.names = F,col.names = F,quote = F)
    k <- length(genes)
    av <- numeric(ncol(expData))
    
    zz <- which(tolower(rownames(expData)) %in% tolower(genes))
    geneExp <- as.matrix(expData[zz, ])
    for (gene in genes) {
      geneExp[gene,] <- (geneExp[gene,] * 10000) / total_counts  
    }
    geneExp <- log2(geneExp+1)
    geneExp <- t(scale(t(geneExp)))
    geneExp[is.nan(geneExp)] <- 0
    av <- av + colSums(geneExp) / k
    
    toPlot <- plotData
    toPlot$Zscore <- NA
    toPlot[colnames(expData), "Zscore"] <- av
    #toPlot$Zscore <- sapply(toPlot$Zscore,function(x){min(x,3)})
    if ("tSNE_1" %in% colnames(plotData)) {
      centers <- toPlot %>% 
        dplyr::group_by_at(opt$column) %>% 
        dplyr::summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) 
      
      pl <- ggplot(data=toPlot[toPlot$Zscore<=quantile(toPlot$Zscore, 0.75),], aes_string(x="tSNE_1", y="tSNE_2")) + theme_classic() + 
        theme(legend.text=element_text(size=5)) +
        expand_limits(x=0, y=0) +
        ggtitle(paste0("Cell type ",ct," - ",sprintf('%.2f%%',avg_props[ct]*100))) + coord_fixed(2/3)
      pl <- pl + geom_point(fill="lightgray", size=0.5, pch=21, stroke=0) +
        geom_point(data=toPlot[toPlot$Zscore>quantile(toPlot$Zscore, 0.75),],
                   aes(fill = Zscore), size=0.75, pch=21, stroke=0) +
        #scale_fill_gradientn(colors=c("lightgray","darkblue","cyan", "lightgreen","yellow",  "red", "red", "darkred","coral4"))
        scale_fill_gradientn(colors=c("lightgray","indianred1","tomato1","red", "red", "darkred")) +
        geom_text(data = centers, aes_string(label=opt$column), size = 3, color = "black")
      
      png(file.path(SAVE.PATH,"SCN",SC_DATASET,"tSNE",paste0("SC_cell_type_",ct,".png")))
      print(pl)
      dev.off()
    }
    
    if ("UMAP_1" %in% colnames(plotData)) {
      centers <- toPlot %>% 
        dplyr::group_by_at(opt$column) %>% 
        dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) 
      
      pl <- ggplot(data=toPlot[toPlot$Zscore<=quantile(toPlot$Zscore, 0.75),], aes_string(x="UMAP_1", y="UMAP_2")) + theme_classic() + 
        theme(legend.text=element_text(size=5)) +
        expand_limits(x=0, y=0) +
        ggtitle(paste0("Cell type ",ct," - ",sprintf('%.2f%%',avg_props[ct]*100))) + coord_fixed(2/3)
      pl <- pl + geom_point(fill="lightgray", size=0.5, pch=21, stroke=0) +
        geom_point(data=toPlot[toPlot$Zscore>quantile(toPlot$Zscore, 0.75),],
                   aes(fill = Zscore), size=0.75, pch=21, stroke=0) +
        #scale_fill_gradientn(colors=c("lightgray","darkblue","cyan", "lightgreen","yellow",  "red", "red", "darkred","coral4"))
        scale_fill_gradientn(colors=c("lightgray","indianred1","tomato1","red", "red", "darkred")) +
        geom_text(data = centers, aes_string(label=opt$column), size = 3, color = "black")
      
      png(file.path(SAVE.PATH,"SCN",SC_DATASET,"UMAP",paste0("SC_cell_type_",ct,".png")))
      print(pl)
      dev.off()
    }
  }
}



rmarkdown::render("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/Docker/docker_linseed_snakemake/app/scripts/reportSCReference.Rmd",
                  output_dir = file.path(opt$save,DT), output_file = paste0(DT,".html"),
                  params = list(sc_dataset = SC_DATASET,
                                res_data = RESULTS_DIR,
                                res_type = "SCN",
                                res_path = file.path(opt$save,DT),
                                min_ct = min_ct,
                                max_ct = max_ct))