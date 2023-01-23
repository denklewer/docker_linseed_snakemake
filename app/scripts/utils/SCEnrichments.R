library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(limma)

set.seed(42)

preparePathwaysGMT <- function(markers,pathway_name,save_path,top_genes=50){
  markers_list <- list()
  for (cluster_ in unique(markers$cluster)) {
    cluster_genes <- (markers %>% filter(cluster == cluster_))$gene
    markers_list[[paste0("sc_",pathway_name,"_",cluster_)]] <- cluster_genes[1:min(top_genes,length(cluster_genes))]
  }
  fgsea::writeGmtPathways(markers_list,gmt.file = save_path)  
  fgsea::gmtPathways(save_path)
}

getSCEnrichments <- function(PATH,DT,sc_pathways,cell_types,
                             calculate_limma=T,filter_=T,filter_limit=0.01) {
  basis <- as.data.frame(fread(file.path(PATH,paste0("ct",cell_types),DT,"best","basis_column.tsv"))) %>% column_to_rownames("V1")
  basis <- basis[-1,]
  if (filter_) {
    props_ <- as.matrix(fread(file.path(PATH,paste0("ct",cell_types),DT,"best","proportions.tsv")) %>% column_to_rownames("V1"))
    avg_props <- apply(props_,1,mean)
    props_ <- props_[avg_props>filter_limit,]
    cell_types <- nrow(props_)
    basis <- basis[,c(paste0("FC_Cell_type_",rownames(props_)),
                      paste0("Cell_type_",rownames(props_)))]
  }
  cols_ <- colnames(basis)[grepl("^Cell.*",colnames(basis))]
  if (calculate_limma) {
    limma_basis <- matrix(0,nrow=nrow(basis),ncol=cell_types)
    colnames(limma_basis) <- paste0("Limma_FC_",cols_)
    for (col_ in cols_) {
      expr <- basis[,cols_]
      condition <- as.integer(col_ != colnames(expr))
      design <- model.matrix(~ condition)
      fit <- lmFit(expr, design)
      fit <- eBayes(fit)
      stats <- topTable(fit, number = nrow(basis))
      limma_basis[,paste0("Limma_FC_",col_)] <- -stats[rownames(basis),"logFC"]
    }
    basis <- cbind(limma_basis,basis)  
  }
  
  
  gseaMatrix <- matrix(0,nrow=length(sc_pathways),ncol=cell_types)
  rownames(gseaMatrix) <- names(sc_pathways)
  colnames(gseaMatrix) <- cols_
  
  fc_column <- "FC_"
  if (calculate_limma) {
    fc_column <- "Limma_FC_"
  }
  for (col_ in cols_) {
    fc_ct <- basis[,paste0(fc_column,col_)]
    names(fc_ct) <- rownames(basis)
    fc_ct_sort <- sort(fc_ct)
    fgseaRes <- fgsea(pathways = sc_pathways, 
                      stats    = fc_ct_sort,
                      scoreType = "pos")
    gseaMatrix[fgseaRes$pathway,col_] <- -log(fgseaRes$pval)*sign(fgseaRes$NES)
  }  
  gseaMatrix
}

# Example of usage
# markers <- as.data.frame(fread("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/FGSEA/pathways/hlca_markers.csv"))
# PATH <- "/Volumes/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/results/GTEX_LUNGS_TPM_HK_CC_F_LIMITX_OMEGA_MAD_1_KNN_Rand_new/09_2022"
# DT <- "20221007_220240"
# cell_types <- 9
# top_genes <- 50
# 
# sc_pathways <- preparePathwaysGMT(markers,"hlca",
#                                   file.path("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/FGSEA/pathways","sc_hlca.gmt"))
# gseaMatrix <- getSCEnrichments(PATH,DT,sc_pathways,cell_types)
# pheatmap(gseaMatrix,scale="none",treeheight_row = F,treeheight_col = F,
#          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))