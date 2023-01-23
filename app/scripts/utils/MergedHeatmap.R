library(ComplexHeatmap)
library(data.table)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(gridExtra)

calculatePairwiseEnrichment <- function(PATH,DT,top_genes=200) {
  cts_ <- sort(as.numeric(gsub("ct","",basename(list.dirs(PATH,recursive = F)))))
  min_ct <- min(cts_)
  max_ct <- max(cts_)
  
  total.merge <- list()
  
  for (ct in min_ct:max_ct) {
    print(ct)
    
    merged.fisher.test <- NULL
    basis_1 <- fread(file.path(PATH,paste0("ct",ct),DT,"best/basis_column.tsv"))
    props_1  <- as.data.frame(basis_1)[1,colnames(basis_1)[grepl("^Cell_",colnames(basis_1))]]
    basis_1 <- as.data.frame(basis_1[-1,]) %>% tibble::column_to_rownames("V1")
    basis_1 <- basis_1[,colnames(basis_1)[grepl("^FC",colnames(basis_1))]]
    
    for (ct_next in min_ct:max_ct) {
      #ct_next <- ct + 1
      basis_2 <- fread(file.path(PATH,paste0("ct",ct_next),DT,"best/basis_column.tsv"))
      props_2  <- as.data.frame(basis_2)[1,colnames(basis_2)[grepl("^Cell_",colnames(basis_2))]]
      basis_2 <- as.data.frame(basis_2[-1,]) %>% tibble::column_to_rownames("V1")
      basis_2 <- basis_2[,colnames(basis_2)[grepl("^FC",colnames(basis_2))]]
      
      test.res <- matrix(0,ct,ct_next)
      rownames(test.res) <- paste0("ct",ct,"_cell_type_",1:ct,"_",props_1)
      colnames(test.res) <- paste0("ct",ct_next,"_cell_type_",1:ct_next,"_",props_2)
      for (ct1 in 1:ct) {
        group_1 <- rownames((basis_1 %>% arrange(desc(.data[[paste0("FC_Cell_type_",ct1)]])))[1:top_genes,])
        for (ct2 in 1:ct_next) {
          group_2 <- rownames((basis_2 %>% arrange(desc(.data[[paste0("FC_Cell_type_",ct2)]])))[1:top_genes,])
          overlap <- intersect(group_1,group_2)
          total <- nrow(basis_1)  
          p.val <- fisher.test(matrix(c(length(overlap), length(group_2)-length(overlap), length(group_1)-length(overlap), total-length(group_2)-length(group_1)+length(overlap)), 2, 2), alternative='greater')$p.value
          test.res[ct1,ct2] <- p.val
        }
      }
      test.res.log <- round(-log10(test.res),3)
      test.res.log[is.infinite(test.res.log)] <- 308.0
      #test.res.log <- as.data.frame(cbind(ct_next,t(test.res.log)))
      test.res.log <- as.data.frame(t(test.res.log))
      test.res.log <- melt(as.matrix(test.res.log))
      test.res.log <- cbind(test.res.log,ct_next)
      colnames(test.res.log)[4] <- c("Var2_ct")
      merged.fisher.test <- rbind(merged.fisher.test,test.res.log)
      
    }
    #total.merge <- cbind(total.merge,merged.fisher.test)
    total.merge[[as.character(ct)]] <- merged.fisher.test
  }
  return(total.merge)
}

filterProportions <- function(data_,filter_limit=0.01){
  selected_rows <- (sapply(rownames(data_),function(x) {
    split_x <- strsplit(x,"_")[[1]]
    as.numeric(split_x[length(split_x)])
  })>filter_limit)
  selected_cols <- (sapply(colnames(data_),function(x) {
    split_x <- strsplit(x,"_")[[1]]
    as.numeric(split_x[length(split_x)])
  })>filter_limit)
  data_[selected_rows,selected_cols]
}

plotHeatmap <- function(data_, from_ct = 0, to_ct = 0,
                        filter_=T, filter_limit=0.01){
  
  if (from_ct == 0) {
    from_ct <- min(as.numeric(names(data_)))
  }
  if (to_ct == 0) {
    to_ct <- max(as.numeric(names(data_)))
  }
  
  full_df <- matrix(0,nrow=0,ncol=3)
  for (ct in from_ct:to_ct) {
    df <- data_[[as.character(ct)]]
    df <- df %>% filter(Var2_ct %in% from_ct:to_ct)
    full_df <- rbind(full_df,df[,c(1,2,3)])
  }
  toPlot <- as.data.frame(dcast(full_df,Var1 ~ Var2))
  toPlot <- toPlot %>% column_to_rownames("Var1")
  if (filter_) {
    toPlot <- filterProportions(toPlot,filter_limit)
  }
  pheatmap(as.matrix(toPlot),scale="none")
}

# Example of use
# PATH <- "/Volumes/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/results/GTEX_MUSCLE_TPM_HK_CC_F_LIMITX_OMEGA_MAD_0_5_KNN_GENES_500_Rand/10_2022"
# DT <- "20221017_230841"
# total.merge <- calculatePairwiseEnrichment(PATH,DT)
# plotHeatmap(total.merge,4,17)




