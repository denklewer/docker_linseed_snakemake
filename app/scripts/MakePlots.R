source("/app/scripts/OutputPlots.R")

metadata_ <- loadRData(snakemake@input[[1]])

png(snakemake@output[[1]])
plotErrors(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()

png(snakemake@output[[2]])
plotNegProportionsChange(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()

png(snakemake@output[[3]])
plotNegBasisChange(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()

png(snakemake@output[[4]])
plotSumToOneChange(metadata_) + labs(title=snakemake@wildcards$"sample")
dev.off()

saveRDS(list(final_Omega=metadata_$final_Omega,
            final_X=metadata_$final_X,
            final_D_w= metadata_$final_D_w,
            final_D_h = metadata_$final_D_h),
                 file=snakemake@output[[5]])
