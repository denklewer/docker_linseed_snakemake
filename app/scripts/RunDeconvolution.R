library(Rcpp)
library(limma)

source('/app/scripts/LinseedMetadata.R')
source('/app/scripts/SinkhornNNLSLinseedC.R')
sourceCpp('/app/scripts/pipeline.cpp')

saveResults <- function(obj_,paths_){
    ## save proportions
      colnames(obj_$full_proportions) <- colnames(obj_$filtered_dataset)
      write.table(obj_$full_proportions,
                  file=paths_[["props_path"]],
                  sep="\t",col.names = NA, row.names = T, quote = F)
    ## save basis row normalized            
      colnames(obj_$full_basis) <- paste0("Cell_type_",1:obj_$cell_types)
      toSave <- obj_$full_basis
      toSave <- obj_$getFoldChange(toSave)
      #toSave <- obj_$getLimmaFoldChange(toSave)
      toSave <- rbind(c(rep(NA,obj_$cell_types),round(apply(obj_$full_proportions,1,mean),4)),toSave)
      rownames(toSave) <- c("avg_proportions",rownames(obj_$filtered_dataset))
      write.table(toSave,file=paths_[["basis_r_path"]],
                  sep="\t",col.names = NA, row.names = T, quote = F)
    ## save basis column normalized
      toSave <- t(t(obj_$full_basis) / (rowSums(t(obj_$full_basis))+1e-10)) 
      toSave <- obj_$getFoldChange(toSave)
      toSave <- toSave[is.nan(toSave)] <- 0
      #toSave <- obj_$getLimmaFoldChange(toSave)
      toSave <- rbind(c(rep(NA,obj_$cell_types),round(apply(obj_$full_proportions,1,mean),4)),toSave)
      rownames(toSave) <- c("avg_proportions",rownames(obj_$filtered_dataset))
      write.table(toSave,file=paths_[["basis_c_path"]],
                  sep="\t",col.names = NA, row.names = T, quote = F)  

    metadata_ <- LinseedMetadata$new(obj_)
    save(metadata_,file=paths_[["meta_path"]])
    
    write.table(obj_$errors_statistics[nrow(obj_$errors_statistics),,drop=FALSE],
                  file=paths_[["stats_path"]],
                  sep="\t",col.names = T, row.names = F, quote = F)  
}

tmp_snk <- readRDS(snakemake@input[["dataset"]])
print(tmp_snk$cell_types)

tmp_snk$readInitValues(snakemake@input[["init_file"]])

blocks_pipeline <- read.csv(snakemake@params[["blocks_pipeline"]],
                            stringsAsFactors=FALSE)

for (idx in 1:nrow(blocks_pipeline)){
    x <- blocks_pipeline[idx,]
    x_limit <- 0
    omega_limit <- 0
    cosine_thresh <- 0
    if ("x_limit" %in% names(x)) {
      x_limit <- as.numeric(x["x_limit"])
    }
    if ("omega_limit" %in% names(x)) {
      omega_limit <- as.numeric(x["omega_limit"])
    }
    if ("cosine_thresh" %in% names(x)) {
      cosine_thresh <- as.numeric(x["cosine_thresh"])
    }
    tmp_snk$runGradientBlock(block_name=x[["block_name"]],
                             coef_der_X = as.numeric(x["coef_der_X"]),
                             coef_der_Omega = as.numeric(x["coef_der_Omega"]),
                             coef_hinge_H = as.numeric(x["coef_hinge_H"]),
                             coef_hinge_W = as.numeric(x["coef_hinge_W"]),
                             coef_pos_D_h = as.numeric(x["coef_pos_D_h"]),
                             coef_pos_D_w = as.numeric(x["coef_pos_D_w"]),
                             iterations = as.numeric(x["iterations"]),
                             limit_X = x_limit,
                             limit_Omega = omega_limit,
                             cosine_thresh = cosine_thresh)
}

saveResults(tmp_snk,list(meta_path=snakemake@output[["meta"]],
                        props_path=snakemake@output[["proportions"]],
                        basis_r_path=snakemake@output[["basis_row"]],
                        basis_c_path=snakemake@output[["basis_column"]],
                        stats_path=snakemake@output[["stats"]]))


