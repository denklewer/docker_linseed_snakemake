library(rhdf5)
library(rjson)
library(jsonlite)
library(Matrix)
library(dplyr)
library(data.table)
library(ggplot2)
library(optparse)

option_list = list(
  make_option(c("--dataset"), type="character", default=NULL),
  make_option(c("--results"), type="character", default=NULL),
  make_option(c("--save"), type="character", default=getwd()),
  make_option(c("--min_ct"), type="numeric", default=NULL),
  make_option(c("--max_ct"), type="numeric", default=NULL),
  make_option(c("--top_genes"), type="numeric", default=100),
  make_option(c("--dt"), type="character", default=NULL),
  make_option(c("--mm"), type="character", default=NULL),
  make_option(c("--column"), type="character", default="Cluster")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


SC_PATH <- opt$dataset #file.path("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/sc_References/SCE/GSE122960_healthy_lung")
SC_DATASET <- basename(SC_PATH) #"GSE131907"
RES_PATH <- opt$results #file.path("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/Snakemake/reports",RESULTS_DIR)
RESULTS_DIR <- basename(RES_PATH)#"LUAD_F_HK_CC_MAD32_GENES_2000"
SAVE.PATH <- opt$save

min_ct <- opt$min_ct
max_ct <- opt$max_ct
top_genes <- opt$top_genes

result <- jsonlite::fromJSON(file.path(SC_PATH,"exp_data.json"))
total_counts <- result$totalCounts
expData <- h5read(file.path(SC_PATH,"data.h5"),"expression")$mat
rownames(expData) <- result$genes
colnames(expData) <- result$barcodes
plotData <- (jsonlite::fromJSON(file.path(SC_PATH,"plot_data.json"))$data) %>% as.data.frame

dir.create(file.path(opt$save,opt$dt),showWarnings = F, recursive = T)

for (cell_types in min_ct:max_ct) {
  
  LOAD.PATH <- file.path(RES_PATH,opt$mm,paste0("ct",cell_types),opt$dt,"best")
  SAVE.PATH <- file.path(opt$save,opt$dt,paste0("ct",cell_types))
  
  basis_ <- data.frame(fread(file.path(LOAD.PATH,"basis_column.tsv"))) %>%
    dplyr::rename(gene = "V1")
  avg_props <- t(basis_[1,-(1:(cell_types+1))])
  basis_ <- basis_[-1,]
  
  dir.create(file.path(SAVE.PATH,"SCE",SC_DATASET,"UMAP"),showWarnings = F, recursive = T)
  dir.create(file.path(SAVE.PATH,"SCE",SC_DATASET,"tSNE"),showWarnings = F, recursive = T)
  
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
    gc()
    
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
      pl <- pl + geom_point(fill="lightgray", size=1.25, pch=21, stroke=0) +
        geom_point(data=toPlot[toPlot$Zscore>quantile(toPlot$Zscore, 0.75),],
                   aes(fill = Zscore), size=1.5, pch=21, stroke=0) +
        scale_fill_gradientn(colors=c("lightgray","indianred1","tomato1","red", "red", "darkred")) +
        geom_text(data = centers, aes_string(label=opt$column), size = 3, color = "black")
        #scale_fill_gradientn(colors=c("darkblue", "cyan", "yellow",  "red", "red", "darkred"))
      
      png(file.path(SAVE.PATH,"SCE",SC_DATASET,"tSNE",paste0("SC_cell_type_",ct,".png")))
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
      pl <- pl + geom_point(fill="lightgray", size=1.25, pch=21, stroke=0) +
        geom_point(data=toPlot[toPlot$Zscore>quantile(toPlot$Zscore, 0.75),],
                   aes(fill = Zscore), size=1.5, pch=21, stroke=0) +
        scale_fill_gradientn(colors=c("lightgray","indianred1","tomato1","red", "red", "darkred")) +
        geom_text(data = centers, aes_string(label=opt$column), size = 3, color = "black")
      
      png(file.path(SAVE.PATH,"SCE",SC_DATASET,"UMAP",paste0("SC_cell_type_",ct,".png")))
      print(pl)
      dev.off()
    }
  }
}



rmarkdown::render("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/Docker/docker_linseed_snakemake/app/scripts/reportSCReference.Rmd",
                  output_dir = file.path(opt$save,opt$dt), output_file = paste0(opt$dt,".html"),
                  params = list(sc_dataset = SC_DATASET,
                                res_data = RESULTS_DIR,
                                res_type = "SCE",
                                res_path = file.path(opt$save,opt$dt),
                                min_ct = min_ct,
                                max_ct = max_ct))