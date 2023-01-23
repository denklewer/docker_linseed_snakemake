library(ggplot2)
library(reshape2)
library(gridExtra)
library(uwot)
library(ggpubr)
library(ggrepel)
library(lsa)

loadRData <- function(fileName){
  #loads an R file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

plotErrors <- function(metadata,variables = c("deconv_error","lamdba_error","beta_error",
    "D_h_error","D_w_error","total_error")) {
      toPlot <- data.frame(metadata$errors_statistics[,variables])
      toPlot$iteration <- 0:(nrow(metadata$errors_statistics)-1)
      toPlot <- melt(toPlot,id.vars="iteration",measure.vars = variables)
      plt <- ggplot(toPlot,aes(x=iteration,y=log10(value),color=variable)) +
      geom_point(size=0.2) +
      geom_line() + theme_minimal()
      plt
}

plotNegBasisChange <- function(metadata_) {
  total_ <- nrow(metadata_$V_row)*metadata_$cell_types
  toPlot <- as.data.frame(metadata_$errors_statistics[,"neg_basis_count",drop=F])
  last_basis_count <- round(toPlot[nrow(toPlot),"neg_basis_count"] / total_,6) * 100
  toPlot$idx <- 1:nrow(toPlot)
  plt <- ggplot(toPlot,aes(y=neg_basis_count,x=idx)) + geom_line() + 
    theme_minimal() + xlab("Iteration") + ylab("Negative basis") +
    annotate("text",  x=Inf, y = Inf, label = paste0(last_basis_count,"%"), vjust=1, hjust=1)
  plt
}

plotNegProportionsChange <- function(metadata_) {
  total_ <- ncol(metadata_$V_row)*metadata_$cell_types
  toPlot <- as.data.frame(metadata_$errors_statistics[,"neg_props_count",drop=F])
  last_prop_count <- round(toPlot[nrow(toPlot),"neg_props_count"] / total_,6) * 100
  toPlot$idx <- 1:nrow(toPlot)
  plt <- ggplot(toPlot,aes(y=neg_props_count,x=idx)) + geom_line() + theme_minimal() + 
    xlab("Iteration") + ylab("Negative proportions")+
    annotate("text",  x=Inf, y = Inf, label = paste0(last_prop_count,"%"), vjust=1, hjust=1)
  plt
}

plotSumToOneChange <- function(metadata_) {
  toPlot <- as.data.frame(metadata_$errors_statistics[,"sum_d_w",drop=F])
  last_sum_dw <- toPlot[nrow(toPlot),"sum_d_w"]
  toPlot$idx <- 1:nrow(toPlot)
  ggplot(toPlot,aes(y=sum_d_w,x=idx)) + geom_line() + theme_minimal() + 
    xlab("Iteration") + ylab("Sum-to-one") +
    annotate("text",  x=Inf, y = Inf, label = sprintf("%.7e",last_sum_dw), vjust=1, hjust=1)
}

plotPoints2D <- function(metadata_,points="init",dims=3) {
      if (!points %in% c("init","current")) {
        stop("Allowed values for points are 'init', 'current'")
      }

      if (points == "init") {
        X <- metadata_$init_X
        Omega <- metadata_$init_Omega
        count_neg_props <- metadata_$init_count_neg_props
        count_neg_basis <- metadata_$init_count_neg_basis
      }
      if (points == "current") {
        X <- metadata_$final_X
        Omega <- metadata_$final_Omega
        count_neg_props <- metadata_$count_neg_props
        count_neg_basis <- metadata_$count_neg_basis
      }

      X <- X[,1:dims]
      Omega <- Omega[1:dims,]

      ## plot X
        toPlot <- as.data.frame(metadata_$V_row %*% t(metadata_$R))[,1:dims]
        colnames(toPlot) <- c("X","Y","Z")
        colnames(X) <- c("X","Y","Z")
        pltX <- ggplot(toPlot, aes(x=Y, y=Z)) +
          geom_point() + 
          geom_polygon(data=as.data.frame(X), fill=NA, color = "green") +
          theme_minimal()
      if (!is.null(count_neg_props)) {
        pltX <- pltX + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_props / (metadata_$cell_types*metadata_$N),4)*100,"%"), vjust=1, hjust=1)
      }

      ## plot Omega  
      toPlot <- as.data.frame(t(metadata_$S %*% metadata_$V_column))[,1:dims]
      colnames(toPlot) <- c("X","Y","Z")
      rownames(Omega) <- c("X","Y","Z")
      pltOmega <- ggplot(toPlot, aes(x=Y, y=Z)) +
        geom_point() + 
        geom_polygon(data=as.data.frame(t(Omega)), fill=NA, color = "green") +
        theme_minimal()
      
      if (!is.null(count_neg_basis)) {
        pltOmega <- pltOmega + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_basis / (metadata_$cell_types*metadata_$M),4)*100,"%"), vjust=1, hjust=1)
      }

      grid.arrange(pltX,pltOmega,nrow=1)
}

plotUMAP <- function(data_, best_run){
  toPlot <- as.data.frame(umap(data_))
  rownames(toPlot) <- rownames(data_)
  toPlot$best <- grepl(paste0("^",best_run),rownames(toPlot))
  toPlot$best_id <- NULL
  init_number <- strsplit(best_run,"_")[[1]]
  init_number <- init_number[length(init_number)]
  toPlot[toPlot$best,"best_id"] <- rownames(toPlot[toPlot$best,])
  toPlot$best_id <- gsub(best_run,init_number,toPlot$best_id)
  ggplot() +  geom_point(data=toPlot[!toPlot$best,],aes(x=V1,y=V2)) + 
    geom_point(data=toPlot[toPlot$best,],aes(x=V1,y=V2,color=best)) +
    geom_text_repel(data=toPlot[toPlot$best,],aes(x=V1,y=V2,label=best_id),
                    min.segment.length = unit(0, 'lines')) +
    theme_minimal() + theme(legend.position = "none") +
    xlab('UMAP1') + ylab('UMAP2')
}

plotCosineUMAP <- function(data_, best_run){
  data_ <- cosine(as.matrix(data_))
  data_[is.nan(data_)] <- 0
  toPlot <- as.data.frame(umap(data_,n_neighbors=2))
  rownames(toPlot) <- rownames(data_)
  toPlot$best <- grepl(paste0("^",best_run),rownames(toPlot))
  toPlot$best_id <- NULL
  init_number <- strsplit(best_run,"_")[[1]]
  init_number <- init_number[length(init_number)]
  toPlot[toPlot$best,"best_id"] <- rownames(toPlot[toPlot$best,])
  toPlot$best_id <- gsub(best_run,init_number,toPlot$best_id)
  ggplot() +  geom_point(data=toPlot[!toPlot$best,],aes(x=V1,y=V2)) + 
    geom_point(data=toPlot[toPlot$best,],aes(x=V1,y=V2,color=best)) +
    geom_text_repel(data=toPlot[toPlot$best,],aes(x=V1,y=V2,label=best_id),
                    min.segment.length = unit(0, 'lines')) +
    theme_minimal() + theme(legend.position = "none") +
    xlab('UMAP1') + ylab('UMAP2')
}

plotAbundance <- function(metadata_){
  colnames(metadata_$full_proportions) <- colnames(metadata_$filtered_dataset)
  rownames(metadata_$full_proportions) <- paste("Cell type",1:metadata_$cell_types)
  toPlot <- melt(metadata_$full_proportions)
  colnames(toPlot) <- c("cell_type","sample","cluster_percent")
  toPlot$cluster_percent <- round(toPlot$cluster_percent*100,2)

  ggbarplot(toPlot, x = "cell_type", y = "cluster_percent",
                fill = "cell_type",  color  ="black",
                position = position_dodge(0.7), add = c("mean_sd", "point"),
                add.params = list(color = "black"), alpha = 0.95) +
  theme(axis.text.x=element_text(angle = 45, hjust=1, size = 12),
        legend.position="none") + xlab("")
}

plotDistribution <- function(data_){
  toPlot <- as.data.frame(c(data_))
  colnames(toPlot) <- "value"
  toPlot$negative <- (toPlot$value < 0)
  ggplot(toPlot) + geom_histogram(aes(x=value,fill=negative),bins = 100) + theme_minimal()

}

plotBasisDistance <- function(data_,basis){
  points <- data_$V_row %*% t(data_$R)
  toPlot <- data.frame(zero_distance=sort(apply(points,1,function(x) sqrt(sum(x^2)))))
  toPlot$cell_type <- "None"
  toPlot$solution <- F
  toPlot$init <- F
  top_genes <- 100
  markers <- list()
  for (ct in 1:data_$cell_types) {
    markers[[paste0("FC_Cell_type_",ct)]] <- c(basis %>% arrange(desc(.data[[paste0("FC_Cell_type_",ct)]])) %>% select(V1))$V1[1:top_genes]
    toPlot[markers[[paste0("FC_Cell_type_",ct)]],"cell_type"] <- paste0("Cell_type_",ct)
  }

  props <- data_$orig_full_proportions / rowSums(data_$orig_full_proportions)
  pointsProp <- props %*% t(data_$R)
  toPlotC <- data.frame(zero_distance=sort(apply(pointsProp,1,function(x) sqrt(sum(x^2)))))
  rownames(toPlotC) <- 1:nrow(toPlotC)
  toPlotC$cell_type <- paste0("Cell_type_",rownames(toPlotC))
  rownames(toPlotC) <- paste0("Solution_",rownames(toPlotC))
  toPlotC$solution <- T
  toPlotC$init <- F
  toPlot <- rbind(toPlot,toPlotC)

  pointsProp <- metadata_$init_X
  toPlotC <- data.frame(zero_distance=sort(apply(pointsProp,1,function(x) sqrt(sum(x^2)))))
  rownames(toPlotC) <- 1:nrow(toPlotC)
  toPlotC$cell_type <- paste0("Cell_type_",rownames(toPlotC))
  rownames(toPlotC) <- paste0("Init_",rownames(toPlotC))
  toPlotC$solution <- F
  toPlotC$init <- T
  toPlot <- rbind(toPlot,toPlotC)

  toPlot <- toPlot[order(toPlot$zero_distance),] 
  toPlot$idx <- 1:nrow(toPlot)
  plt <- ggplot(toPlot) + geom_point(aes(x=idx,y=zero_distance,color=cell_type)) + 
  geom_point(data=toPlot %>% filter(solution),aes(x=idx,y=zero_distance),color="black") +
  geom_point(data=toPlot %>% filter(init),aes(x=idx,y=zero_distance),color="blue") +
  geom_vline(data=toPlot %>% filter(solution),aes(xintercept=idx),color="black") + 
  geom_vline(data=toPlot %>% filter(init),aes(xintercept=idx),color="blue") + 
  theme_minimal() + facet_wrap(~cell_type)
  plt
}

plotProportionsDistance <- function(data_,proportions) {
  points <- t(data_$S %*% data_$V_column)

  toPlot <- data.frame(zero_distance=sort(apply(points,1,function(x) sqrt(sum(x^2)))))
  toPlot$cell_type <- "None"
  toPlot$solution <- F

  basis_ <- t(t(data_$orig_full_basis) / rowSums(t(data_$orig_full_basis)))
  pointsBasis <- t(data_$S %*% basis_)
  toPlotC <- data.frame(zero_distance=sort(apply(pointsBasis,1,function(x) sqrt(sum(x^2)))))
  rownames(toPlotC) <- 1:nrow(toPlotC)
  toPlotC$cell_type <- paste0("Cell_type_",rownames(toPlotC))
  rownames(toPlotC) <- paste0("Solution_",rownames(toPlotC))
  toPlotC$solution <- T
  toPlot <- rbind(toPlot,toPlotC)

  toPlot <- toPlot[order(toPlot$zero_distance),] 
  toPlot$idx <- 1:nrow(toPlot)

  plt <- ggplot(toPlot) + geom_point(aes(x=idx,y=zero_distance),color="lightgray") + 
    geom_point(data=toPlot %>% filter(solution),aes(x=idx,y=zero_distance,color=cell_type)) +
    geom_text_repel(data=toPlot %>% filter(solution),aes(x=idx,y=zero_distance,label=cell_type),
                    max.overlaps=100) +
    theme_minimal()
  plt
}

plotCosineHeatmap <- function(data_){
  
  rownames(data_) <- as.character(1:nrow(data_))
  colnames(data_) <- as.character(1:ncol(data_))
  toPlot <- round(cosine(data_),2)
  toPlot[is.nan(toPlot)] <- 0
  pheatmap::pheatmap(toPlot,scale="none",cluster_cols = F, cluster_rows = F, display_numbers = T,
                    show_rownames = T,show_colnames = T)
}