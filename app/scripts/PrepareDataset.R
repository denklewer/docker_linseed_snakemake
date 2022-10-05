library(Rcpp)
source("/app/scripts/SinkhornNNLSLinseedC.R")
source("/app/scripts/PreprocessingPlots.R")
sourceCpp('/app/scripts/pipeline.cpp')

filterZeroMAD <- function(dataset){
    mad_ <- apply(dataset,1,mad)
    dataset[mad_>0,]
}

data_ <- readRDS(snakemake@params[["dataset"]])
data_ <- filterZeroMAD(data_)
tmp_snk <- SinkhornNNLSLinseed$new(dataset = snakemake@config[["dataset"]], 
                                    path = snakemake@output[[1]], 
                                    data = data_,
                                    analysis_name = snakemake@config[["analysis_name"]],
                                    cell_types = snakemake@config[["cell_types"]],
                                    k_genes = snakemake@config[["k_genes"]],
                                    k_samples = snakemake@config[["k_samples"]])
print(tmp_snk$cell_types)
if (!is.null(snakemake@config[["top_median"]])){
    tmp_snk$metric <- "median"
    med_limit <- min(tmp_snk$M,snakemake@config[["top_median"]])
    top_genes <- names(sort(tmp_snk$genes_median,decreasing = T)[1:med_limit])
    png(snakemake@output[["top_mad"]])
    print(plotTopMADGenes(tmp_snk,top_genes))
    dev.off()
    tmp_snk$selectTopGenes(med_limit)
} else if (!is.null(snakemake@config[["top_mad"]])){
    tmp_snk$metric <- "mad"
    mad_limit <- min(tmp_snk$M,snakemake@config[["top_mad"]])
    top_genes <- names(sort(tmp_snk$genes_mad,decreasing = T)[1:mad_limit])
    png(snakemake@output[["top_mad"]])
    print(plotTopMADGenes(tmp_snk,top_genes))
    dev.off()
    tmp_snk$selectTopGenes(mad_limit)
} else {
    min_mad = 0
    min_median = 0
    max_mad = Inf
    max_median = Inf
    if (!is.null(snakemake@config[["min_mad"]])){
        min_mad=snakemake@config[["min_mad"]]
    }
    if (!is.null(snakemake@config[["max_mad"]])){
        max_mad=snakemake@config[["max_mad"]]
    }
    if (!is.null(snakemake@config[["min_median"]])){
        min_median=snakemake@config[["min_median"]]
    }
    if (!is.null(snakemake@config[["max_median"]])){
        max_median=snakemake@config[["max_median"]]
    }
    print(paste0(min_mad,"-",max_mad))
    print(paste0(min_median,"-",max_median))
    top_genes1 <- names(tmp_snk$genes_mad[(tmp_snk$genes_mad >= min_mad) & (tmp_snk$genes_mad <= max_mad)])
    top_genes2 <- names(tmp_snk$genes_median[(tmp_snk$genes_median >= min_median) & (tmp_snk$genes_median <= max_median)])
    top_genes <- intersect(top_genes1,top_genes2)
    png(snakemake@output[["top_mad"]])
    print(plotTopMADGenes(tmp_snk,top_genes1))
    dev.off()

    png(snakemake@params[["top_mad_joined"]])
    print(plotTopMADGenes(tmp_snk,top_genes))
    dev.off()

    png(snakemake@params[["mad_med_scatter"]])
    print(plotMADExpMedianScatter(tmp_snk,top_genes))
    dev.off()
    
    tmp_snk$filterByMADExpMedian(min_mad,max_mad,min_median,max_median)
}
cat(paste0("Genes: ",tmp_snk$M),paste0("Samples: ",tmp_snk$N),
        file=snakemake@output[["metadata"]],sep="\n")
tmp_snk$scaleDataset(snakemake@config[["scale_iterations"]])
png(snakemake@output[["svd_before"]])
plotSVD(tmp_snk,snakemake@output[["svd_before_plot"]])
dev.off()
tmp_snk$getSvdProjectionsNew()
tmp_snk$calculateDistances()

png(snakemake@output[["knn_distance_before"]])
plotKNNDistances(tmp_snk,snakemake@config[["thresh_genes"]],snakemake@config[["thresh_samples"]])
dev.off()

tmp_snk$filterByKNN(thresh_genes = snakemake@config[["thresh_genes"]],
                    thresh_samples = snakemake@config[["thresh_samples"]],
                    iterations = snakemake@config[["scale_iterations"]])

tmp_snk$calculateDistances()
png(snakemake@output[["knn_distance_after"]])
plotKNNDistances(tmp_snk,snakemake@config[["thresh_genes"]],snakemake@config[["thresh_samples"]])
dev.off()

png(snakemake@output[["distance_before"]])
plotDistances(tmp_snk,snakemake@config[["filter_genes"]],snakemake@config[["filter_samples"]])
dev.off()

tmp_snk$filterByDistance(filter_genes = snakemake@config[["filter_genes"]],
                         filter_samples = snakemake@config[["filter_samples"]],
                         iterations = snakemake@config[["scale_iterations"]])
tmp_snk$calculateDistances()                         
png(snakemake@output[["distance_after"]])
plotDistances(tmp_snk,0,0)
dev.off()
png(snakemake@output[["svd_after"]])
plotSVDMerged(tmp_snk,snakemake@output[["svd_before_plot"]])
dev.off()
saveRDS(tmp_snk,file=snakemake@output[[1]])