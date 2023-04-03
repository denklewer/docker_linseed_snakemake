library(R6)
library(linseed)
library(NMF)
library(ggplot2)
library(combinat)
library(progress)
library(corpcor)
library(MASS)
library(nnls)
library(gridExtra)
library(dbscan)

SinkhornNNLSLinseed <- R6Class(
  "SinkhornNNLSLinseed",
  public = list(
    filtered_samples = NULL,
    filtered_dataset = NULL,
    raw_dataset = NULL,
    linseed_object = NULL,
    path_ = NULL,
    analysis_name = NULL,
    dataset = NULL,
    topGenes = NULL,
    cell_types = NULL,
    samples = NULL,
    coef_der_X = NULL,
    coef_der_Omega = NULL,
    coef_hinge_H = NULL,
    coef_hinge_W = NULL,
    coef_pos_D_w = NULL,
    coef_pos_D_h = NULL,
    global_iterations = NULL,
    new_points = NULL,
    new_samples_points = NULL,
    restore_all_X = NULL,
    restore_all_Omega = NULL,
    orig_full_proportions = NULL,
    full_proportions = NULL,
    orig_full_basis = NULL,
    full_basis = NULL,
    
    preprocessing_cell_types = NULL,
    preprocessing_organism = NULL,
    preprocessing_coding_genes = NULL,
    
    plane_distance_genes = NULL,
    plane_distance_samples = NULL,
    zero_distance_genes = NULL,
    zero_distance_samples = NULL,
    
    mean_radius_X = NULL,
    mean_radius_Omega = NULL,
    filters_pipeline = matrix(0,nrow=0,ncol=3),
    data = NULL,
    V_row = NULL,
    V_column = NULL,
    M = NULL,
    N = NULL,
    Sigma = NULL,
    R = NULL,
    S = NULL,
    R_ext = NULL,
    S_ext = NULL,
    A = NULL,
    B = NULL,
    X = NULL,
    Omega = NULL,
    W_ = NULL,
    H_ = NULL,

    D_h = NULL,
    D_w = NULL,
    D = NULL,
    D_v_row = NULL,
    D_v_column = NULL,

    init_D_w = NULL,
    init_D_h = NULL,
    init_D = NULL,
    init_X = NULL,
    init_H = NULL,
    init_W = NULL,
    init_Omega = NULL,
    
    init_count_neg_props = NULL,
    init_count_neg_basis = NULL,
    count_neg_props = NULL,
    count_neg_basis = NULL,
    errors_statistics = NULL,
    points_statistics_X = NULL,
    points_statistics_Omega = NULL,
    blocks_statistics = NULL,
    init_errors_statistics = NULL,
    genes_mean = NULL,
    genes_sd = NULL,
    genes_mad = NULL,
    genes_median = NULL,
    top_genes = NULL,
    metric = NULL,
    
    linearizeDataset = function(ge) {
      if (self$is_logscale(ge))
        return(2^ge - 1)
      return(ge)
    },
    
    logDataset = function(ge) {
      if (self$is_logscale(ge))
        return(ge)
      return(log2(ge + 1))
    },
    
    is_logscale = function(x) {
      qx <- quantile(as.numeric(x), na.rm = T)
      if (qx[5] - qx[1] > 100 || qx[5] > 100) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    },
    
    getFoldChange = function(signatures) {
      cell_types_fc <-
        matrix(0,
               ncol = ncol(signatures),
               nrow = nrow(signatures))
      rownames(cell_types_fc) <-
        rownames(signatures)
      colnames(cell_types_fc) <-
        paste0("FC_", colnames(signatures))
      for (i in 1:ncol(cell_types_fc)) {
        cell_types_fc[, i] <- apply(signatures, 1, function(x) {
          x[i] / mean(x[-i])
        })
      }
      cbind(cell_types_fc, signatures)
    },
    
    getLimmaFoldChange = function(signatures){
      limma_fc <- matrix(0,nrow=nrow(signatures),ncol=self$cell_types)
      colnames(limma_fc) <- paste0("LimmaFC_", 1:self$cell_types)
      for (ct in 1:self$cell_types) {
        expr <- signatures[,grepl("^Cell.*",colnames(signatures))]
        condition <- as.integer(paste0("Cell_type_",ct) != colnames(expr))
        design <- model.matrix(~ condition)
        fit <- lmFit(expr, design)
        fit <- eBayes(fit)
        stats <- topTable(fit, number = nrow(signatures))
        limma_fc[,ct] <- -stats[rownames(signatures),"logFC"]
      }
      cbind(limma_fc, signatures)
    },
    
    fcnnls_coefs = function(H_0, filtered_dataset) {
      W <- (.fcnnls(H_0, filtered_dataset, pseudo = TRUE))$coef
      W
    },
    
    preprocessing = function(dataset) {
      # Filter functions applied to rows (gene names), return new genes or boolean list
      conditions  <- c(
        ## remove rows with NaN values
        "NaN containing" = function(dataset) complete.cases(dataset),
        ## filter out non-coding genes
        "non-coding" = function(dataset) {
            if (is.null(self$preprocessing_coding_genes)) {
                if (self$preprocessing_organism == "Mouse") {
                gene_names <- readRDS("/app/scripts/mice_coding_genes_v2.rds")
                }
                else {
                gene_names <- readRDS("/app/scripts/coding_genes_v2.rds")
                }
            }
            else {
                gene_names <- self$preprocessing_coding_genes
            }

            intersect(rownames(dataset), gene_names)
          },
        ## filter out RPL/RPS genes
        "RPL/RPS" =  function(dataset) !grepl("^(RPL|RPS).+", rownames(dataset)),
        ## filter out LOC genes
        "LOC" = function(dataset) !grepl("^^LOC\\d+", rownames(dataset)),
        ## filter out C...orf genes
        "C...orf" = function(dataset)!grepl("^C\\w+orf\\d+", rownames(dataset)),
        ## filter out zero mad genes
        "zero mad" = function(dataset) {
          mad_genes <- apply(dataset, 1, mad)
          names(mad_genes[mad_genes > 0])
          }
      )
      samples <- ncol(dataset)
      count_before <- nrow(dataset)
      print(paste("Genes before filtering:", count_before, "genes"))
      self$filters_pipeline <- rbind(self$filters_pipeline, c("Original data", count_before, samples))

      # Apply all filters sequentially
      for (condition_name in names(conditions)) {
        condition <- conditions[[condition_name]]
        count_before <- nrow(dataset)
        dataset <- dataset[condition(dataset), ]
        count_after <- nrow(dataset)
        filtered_count <- count_before - count_after
        print(paste("Removed:", filtered_count, condition_name , "genes"))
        self$filters_pipeline <- rbind(self$filters_pipeline, c(paste("After filtering",condition_name, "genes"), count_after, samples))
      }
      count_after <- nrow(dataset)
      print(paste("Genes after preprocessing:", count_after, "genes left"))
      dataset
    },

    initialize = function(dataset,
                          path,
                          analysis_name,
                          cell_types,
                          filtered_samples = c(),
                          topGenes = 100000,
                          data = NULL,
                          metric="mad",
                          coef_der_X = 0.001,
                          coef_der_Omega = 0.001,
                          coef_pos_D_w = 0.01,
                          coef_pos_D_h = 0.01,
                          coef_hinge_H = 100,
                          coef_hinge_W = 10,
                          global_iterations = 10000,
                          data_type="bulk",
                          annotate = F,
                          geneSymbol = "Gene symbol",
                          linearize = F,
                          preprocessing = F,
                          preprocessing_cell_types=20,
                          preprocessing_organism = "Human",
                          preprocessing_coding_genes = NULL
                          ) {
      self$filtered_samples <- filtered_samples
      self$dataset <- dataset
      self$path_ <- path
      self$analysis_name <- analysis_name
      self$topGenes <- topGenes
      self$cell_types <- cell_types
      self$preprocessing_cell_types <- preprocessing_cell_types
      self$preprocessing_organism <- preprocessing_organism
      self$preprocessing_coding_genes <- preprocessing_coding_genes
      
      self$data <- data
      
      self$coef_der_X <- coef_der_X
      self$coef_der_Omega <- coef_der_Omega
      self$coef_hinge_H <- coef_hinge_H
      self$coef_hinge_W <- coef_hinge_W
      self$coef_pos_D_w <- coef_pos_D_w
      self$coef_pos_D_h <- coef_pos_D_h
      
      self$global_iterations <- global_iterations
      
      if (!is.null(data)) {
        input_data <- self$data
        
        if (annotate) {
          fdata <- input_data[, geneSymbol, drop = FALSE]
          input_data <- input_data[,-which(colnames(input_data)==geneSymbol)]
        }
      } else {
        gse <- getGEO(self$dataset, AnnotGPL = T)
        if (length(gse) > 1) {
          stop("This GSE has multiple expression sets. It's probably multiseries. Provide single series experiment")
        }
        gse <- gse[[1]]
        input_data <- Biobase::exprs(gse)
        if (annotate) {
          fdata <- fData(gse)[, geneSymbol, drop = FALSE]
        }
        
      }
      
      if (annotate) {
        probes <- setNames(as.character(fdata[, 1]), rownames(fdata))
        probes <- probes[!(is.na(probes) | probes=="")]

        geneToProbes <- lapply(split(probes, probes), names)
        
        bestProbes <- sapply(geneToProbes, function(probes) {
          subset <- input_data[probes, , drop = FALSE]
          probes[which.max(rowMeans(subset))]
        })
        input_data <- input_data[bestProbes, ]
        rownames(input_data) <- names(geneToProbes)
      }
      
      if (length(self$filtered_samples) != 0) {
        input_data <- input_data[,self$filtered_samples]
      }
      
      if (preprocessing) {
        input_data <- self$preprocessing(input_data)  
      }
      
      if (data_type == "microarray") {
        input_data <- self$logDataset(input_data)
        input_dataCopy <- normalize.quantiles(input_data)
        colnames(input_dataCopy) <- colnames(input_data)
        rownames(input_dataCopy) <- rownames(input_data)
        input_data <- input_dataCopy
      }
      
      if (linearize) {
        input_data <- self$linearizeDataset(input_data)  
      }
      
      self$raw_dataset <- as.matrix(input_data)
      self$filtered_dataset <- self$raw_dataset / rowSums(self$raw_dataset)
      
      self$genes_mean <- apply(self$raw_dataset,1,mean)
      self$genes_median <- apply(self$raw_dataset,1,median)
      self$genes_sd <- apply(self$raw_dataset,1,sd)
      self$genes_mad <- apply(self$raw_dataset,1,mad)
      
      self$samples <- ncol(self$filtered_dataset)
      self$metric <- metric

      self$N <- ncol(self$filtered_dataset)
      self$M <- nrow(self$filtered_dataset)

      self$filters_pipeline <- rbind(self$filters_pipeline,c("Algorithm input data",self$M,self$N))
      
    },
    
    filterByMAD = function(min_mad=-Inf,max_mad=Inf){
      self$top_genes <- names(self$genes_mad[(self$genes_mad >= min_mad) & (self$genes_mad <= max_mad)])
      self$filtered_dataset <- self$filtered_dataset[self$top_genes,]
      self$M <- nrow(self$filtered_dataset)
    },

    filterByMADExpMedian = function(min_mad=-Inf,max_mad=Inf,min_median=-Inf,max_median=Inf){
      top_genes1 <- names(self$genes_mad[(self$genes_mad >= min_mad) & (self$genes_mad <= max_mad)])
      top_genes2 <- names(self$genes_median[(self$genes_median >= min_median) & (self$genes_median <= max_median)])
      self$top_genes <- intersect(top_genes1,top_genes2)
      self$filtered_dataset <- self$filtered_dataset[self$top_genes,]
      self$M <- nrow(self$filtered_dataset)
    },
    
    selectTopGenes = function(genes_number=10000) {
      if (self$metric == "mean") {
        dataset_ <- self$genes_mean
      } else if (self$metric == "sd") {
        dataset_ <- self$genes_sd
      } else if (self$metric == "mad") {
        dataset_ <- self$genes_mad
      } else if (self$metric == "median") {
        dataset_ <- self$genes_median
      } else {
        stop("Metric not find. Available metrics: mean, sd, mad")
      }
      self$top_genes <- names(sort(dataset_,decreasing = T)[1:genes_number])
      self$filtered_dataset <- self$filtered_dataset[self$top_genes,]
      self$M <- nrow(self$filtered_dataset)
    },
    
    scaleDataset = function(iterations = 20){
      V <- self$raw_dataset[rownames(self$filtered_dataset),
                            colnames(self$filtered_dataset)]
      
      scaled <- scaleDataset(V,iterations)
      self$V_row <- scaled$V_row
      rownames(self$V_row) <- rownames(self$filtered_dataset)
      colnames(self$V_row) <- colnames(self$filtered_dataset)
      
      self$V_column <- scaled$V_column
      rownames(self$V_column) <- rownames(self$filtered_dataset)
      colnames(self$V_column) <- colnames(self$filtered_dataset)
    },
    
    getSvdProjectionsNew = function(k = self$cell_types){
      svd_ <- svd(self$V_row)
      self$S <- t(svd_$u[,1:k])
      self$R <- t(svd_$v[,1:k])
      self$Sigma <- diag(svd_$d[1:k])
      if (all(self$R[1,]<0)) {
        self$S[1,] <- -self$S[1,]  
        self$R[1,] <- -self$R[1,]
      }
      self$A <- matrix(apply(self$R,1,sum),ncol=1,nrow=self$cell_types)
      self$new_points <- self$V_row %*% t(self$R)
      
      self$B <- matrix(apply(self$S,1,sum),ncol=1,nrow=self$cell_types)
      self$new_samples_points <- t(self$S %*% self$V_column)
      
      self$mean_radius_X <- mean(apply(self$new_points[,-1],1,function(x){norm(x,"2")}))
      self$mean_radius_Omega <- mean(apply(self$new_samples_points[,-1],1,function(x){norm(x,"2")}))
      
      dims <- 1:min(nrow(self$V_row), ncol(self$V_row))
      self$S_ext <- t(svd_$u[, dims])
      self$R_ext <- t(svd_$v[, dims])
      if (all(self$R_ext[1,]<0)) {
        self$S_ext[1,] <- -self$S_ext[1,]
        self$R_ext[1,] <- -self$R_ext[1,]
      }
      
      self$restore_all_X <- self$V_row %*% t(self$R_ext)
      self$restore_all_Omega <- t(self$S_ext %*% self$V_column)
    },
    
    calculatePartialDistances = function(data, with_dims) {
      data_flt <- data
      data_flt[, -with_dims] <- 0
      sqrt(rowSums(data_flt^2))
    },
    
    calculateDistances = function() {
      if (is.null(self$R)) {
        stop("Run getSvdProjectionsNew first")
      }

      self$plane_distance_genes <- sort(self$calculatePartialDistances(self$restore_all_X, 
                                                                  with_dims = -c(1, 2:self$preprocessing_cell_types)),decreasing=T)
      self$plane_distance_samples <- sort(self$calculatePartialDistances(self$restore_all_Omega, 
                                                                    with_dims = -c(1, 2:self$preprocessing_cell_types)),decreasing=T)

      self$zero_distance_genes <- sort(self$calculatePartialDistances(self$restore_all_X, 
                                                                      with_dims = 2:self$preprocessing_cell_types))
      self$zero_distance_samples <- sort(self$calculatePartialDistances(self$restore_all_Omega, 
                                                                        with_dims = 2:self$preprocessing_cell_types))
    },
    
    distanceCutoff = function(distances, threshold, samples = T) {
      keep_names <- names(distances[distances < threshold])
      if (samples) {
        self$filtered_dataset <- self$filtered_dataset[,keep_names]
        self$N <- ncol(self$filtered_dataset)
      } else {
        self$filtered_dataset <- self$filtered_dataset[keep_names,]
        self$M <- nrow(self$filtered_dataset)
      }
    },
    
    thresholdCutoff = function(
      method = "n_sigma",
      distance = "zero_distance",
      samples = T,
      n_sigma = 3,
      zd_quantile = 0.9,
      recalculate_distances = TRUE
    ) {
      if (samples) {
        metric <- paste0(distance,"_samples")
      } else {
        metric <- paste0(distance,"_genes")
      }
      distances <- sort(get(metric,self))
      if (method == "n_sigma") {
        middle <- mean(distances)
        upper_bound <- middle + n_sigma * sd(distances)
      } else if (method == "zd_quantile") {
        upper_bound <- quantile(distances, zd_quantile)
      }
      self$distanceCutoff(distances, upper_bound, samples)
      if (recalculate_distances) {
        self$scaleDataset()
        self$getSvdProjectionsNew()
        self$calculateDistances()
      }
    },
    
    removeOutliers = function(cutoff_samples = T, cutoff_genes = T, filter_by_plane = T){
      if (cutoff_samples + cutoff_genes > 0) {
        prev_length_samples <- Inf
        prev_length_genes <- Inf
        while(cutoff_samples && (prev_length_samples > self$N) || cutoff_genes && (prev_length_genes > self$M)) {
          if (cutoff_samples) {
            prev_length_samples <- self$N
              new_anno <- self$thresholdCutoff(method = "n_sigma", samples = T, recalculate_distances = !cutoff_genes)
            }
            if (cutoff_genes) {
              prev_length_genes <- self$M
              new_anno <- self$thresholdCutoff(method = "n_sigma", samples = F)
            }
          }
      }
      if (filter_by_plane) {
        if (cutoff_samples) {
          self$thresholdCutoff(method = "zd_quantile", distance = "plane_distance", samples = F)
        }
        if (cutoff_genes) {
          self$thresholdCutoff(method = "zd_quantile", distance = "plane_distance", samples = T)
          }
      }
      },
    
    
    
    ## Initializations
    
    selectInitOmega = function(seed = NULL,
                               points = self$new_samples_points) {
      set.seed(seed)
      
      restored <- t(self$S) %*% t(points)
      p <- self$cell_types
      x <- t(points)
      u <- rowMeans(x)
      y <- x / matrix(kronecker(colSums(x * u), rep(1, p)), nrow=p)
      
      indice <- rep(0, p)
      A <- matrix(0, nrow=p, ncol=p)
      A[p, 1] <- 1
      for (i in 1:p) {
        w <- matrix(runif(p), ncol=1)
        f <- w - A %*% pseudoinverse(A) %*% w
        f <- f / sqrt(sum(f^2))

        v <- t(f) %*% y
        indice[i] <- which.max(abs(v))
        A[, i] <- y[, indice[i]]
      }
      Ae <- restored[, indice]
      ## Omega
      self$init_Omega <- self$S %*% Ae

      ## D
      self$init_D_w <- ginv(self$init_Omega) %*% self$B
      self$init_D_h <- self$init_D_w * (self$N/self$M)
      
      ## X
      V__ <- self$S %*% self$V_row %*% t(self$R)
      self$init_X <- ginv(self$init_Omega %*% diag(self$init_D_w[,1])) %*% V__
    },

    selectInitX = function(seed = NULL,
                           points = self$new_points) {
      set.seed(seed)
      
      restored <- points %*% self$R
      p <- self$cell_types

      x <- t(points)
      u <- rowMeans(x)
      y <- x / matrix(kronecker(colSums(x * u), rep(1, p)), nrow=p)

      indice <- rep(0, p)
      A <- matrix(0, nrow=p, ncol=p)
      A[p, 1] <- 1
      for (i in 1:p) {
        w <- matrix(runif(p), ncol=1)
        f <- w - A %*% pseudoinverse(A) %*% w;
        f <- f / sqrt(sum(f^2))

        v <- t(f) %*% y
        indice[i] <- which.max(abs(v))
        A[, i] <- y[, indice[i]]
      }
      Ae <- restored[indice, ]
      ## X
      self$init_X <- Ae %*% t(self$R)
      ## D
      self$init_D_h <- ginv(t(self$init_X)) %*% self$A
      self$init_D_w <- self$init_D_h * (self$M/self$N)
      ## Omega
      V__ <- self$S %*% self$V_row %*% t(self$R)
      self$init_Omega <- V__ %*% ginv(diag(self$init_D_w[,1]) %*% self$init_X)

    },
    
    selectInitXConvex = function(r_tilda=0.95){
      limit_num_ <- floor(nrow(self$V_row)*r_tilda)
      self$init_X <- diag(apply(self$new_points[,-1],2,function(x){
        sort(x)[limit_num_]
      }))
      self$init_X <- cbind(1/sqrt(self$N),
                              rbind(self$init_X,
                                    -1/(self$cell_types-1) * apply(self$init_X,1,sum)))
      self$init_D_h <- ginv(t(self$init_X)) %*% self$A
      self$init_D_w <- self$init_D_h * (self$M/self$N)
      ## Omega
      V__ <- self$S %*% self$V_row %*% t(self$R)
      self$init_Omega <- V__ %*% ginv(diag(self$init_D_w[,1]) %*% self$init_X)
    },
    
    selectInitOmegaConvex = function(r_tilda=0.95){
      limit_num_ <- floor(ncol(self$V_row)*r_tilda)
      self$init_Omega <- diag(apply(self$new_samples_points[,-1],2,function(x){
        sort(x)[limit_num_]
      }))
      self$init_Omega <- t(cbind(1/sqrt(self$M),
                                    rbind(self$init_Omega,
                                          -1/(self$cell_types-1) * apply(self$init_Omega,1,sum))))
      ## D
      self$init_D_w <- ginv(self$init_Omega) %*% self$B
      self$init_D_h <- self$init_D_w * (self$N/self$M)
      ## X
      V__ <- self$S %*% self$V_row %*% t(self$R)
      self$init_X <- ginv(self$init_Omega %*% diag(self$init_D_w[,1])) %*% V__
    },
    
    selectInitXSubset = function(thresh=2000) {
      N <- max(length(self$zero_distance_genes),thresh)
      genes_subset <- names(self$zero_distance_genes[(N-thresh):N])
      points_subset <- self$new_points[genes_subset,]
      print(dim(points_subset))
      self$selectInitX(points=points_subset)
    },

    selectInitOmegaSubset = function(thresh=100) {
      M <- max(length(self$zero_distance_samples),thresh)
      samples_subset <- names(self$zero_distance_samples[(M-thresh+1):M])
      points_subset <- self$new_samples_points[samples_subset,]
      self$selectInitOmega(points=points_subset)
    },
    
    selectInitRandom = function(seed = NULL,
                                n = 1000) {
      set.seed(seed)
      
      idxTableX <- matrix(0,ncol=self$cell_types+1,nrow=n)
      idxTableOmega <- matrix(0,ncol=self$cell_types+1,nrow=n)
      for (i in 1:n) {
        #Omega
        ids_Omega <- sample(1:self$N,self$cell_types)  
        Ae <- self$V_column[,ids_Omega]
        init_Omega <- self$S %*% Ae
        metric_Omega <- sqrt(sum(apply(init_Omega[-1,],1,mean)^2))
        
        idxTableOmega[i,] <- c(ids_Omega,metric_Omega)
        
        
        #X
        ids_X <- sample(1:self$M,self$cell_types)
        Ae <- self$V_row[ids_X,]
        init_X <- Ae %*% t(self$R)
        metric_X <- sqrt(sum(apply(init_X[,-1],2,mean)^2))
        
        idxTableX[i,] <- c(ids_X,metric_X)
      }
      
      idxTableOmega <- idxTableOmega[order(idxTableOmega[,(self$cell_types+1)],decreasing=F),]
      idxTableX <- idxTableX[order(idxTableX[,(self$cell_types+1)],decreasing=F),]
      
      #Omega
      ids_Omega <- idxTableOmega[1,,drop=F]
      Ae <- self$V_column[,ids_Omega]
      self$init_Omega <- self$S %*% Ae
      
      #X
      ids_X <- idxTableX[1,,drop=F]
      Ae <- self$V_row[ids_X,]
      self$init_X <- Ae %*% t(self$R)
      
      V__ <- self$S %*% self$V_row %*% t(self$R)
      ## calculate D_w and D_h
      ## vectorizing deconvolution
      vec_mtx <- matrix(0,self$cell_types*self$cell_types,self$cell_types)
      for (col_ in 1:self$cell_types) {
        vec_mtx[,col_] <- cbind(c(t(t(self$init_Omega[,col_])) %*% self$init_X[col_,]))
      }
      ## adding sum-to-one constraint
      self$init_D_w <- matrix(nnls(rbind(vec_mtx,self$init_Omega),rbind(cbind(c(V__)),self$B))$x,nrow=self$cell_types,ncol=1)
      self$init_D_h <- self$init_D_w * (self$N/self$M)
    },
    
    selectInitRandomCentered = function(seed = NULL) {
      set.seed(seed)
      
      #Omega
      ids_Omega <- sample(1:self$N,(self$cell_types-1))  
      Ae <- self$V_column[,ids_Omega]
      
      init_Omega <- self$S %*% Ae
      self$init_Omega <- cbind(c(1/sqrt(self$M),-apply(init_Omega[-1,],1,sum)),
                          init_Omega)
      
      #X
      ids_X <- sample(1:self$N,(self$cell_types-1))
      Ae <- self$V_row[ids_X,]
      init_X <- Ae %*% t(self$R)
      self$init_X <- rbind(c(1/sqrt(self$N),-apply(init_X[,-1],2,sum)),
                           init_X)
      
      V__ <- self$S %*% self$V_row %*% t(self$R)
      ## calculate D_w and D_h
      ## vectorizing deconvolution
      vec_mtx <- matrix(0,self$cell_types*self$cell_types,self$cell_types)
      for (col_ in 1:self$cell_types) {
        vec_mtx[,col_] <- cbind(c(t(t(self$init_Omega[,col_])) %*% self$init_X[col_,]))
      }
      ## adding sum-to-one constraint
      self$init_D_w <- matrix(nnls(rbind(vec_mtx,self$init_Omega),rbind(cbind(c(V__)),self$B))$x,nrow=self$cell_types,ncol=1)
      self$init_D_h <- self$init_D_w * (self$N/self$M)
    },

  initWithSubset = function(n,top) {
      idxTableX <- matrix(0,ncol=self$cell_types+1,nrow=n)
      idxTableOmega <- matrix(0,ncol=self$cell_types+1,nrow=n)
      for (i in 1:n) {
        #Omega
        ids_Omega <- sample(1:self$N,self$cell_types)  
        Ae <- self$V_column[,ids_Omega]
        init_Omega <- self$S %*% Ae
        metric_Omega <- sqrt(sum(apply(init_Omega[-1,],1,mean)^2))
        
        idxTableOmega[i,] <- c(ids_Omega,metric_Omega)
        
        
        #X
        ids_X <- sample(1:self$N,self$cell_types)
        Ae <- self$V_row[ids_X,]
        init_X <- Ae %*% t(self$R)
        metric_X <- sqrt(sum(apply(init_X[,-1],2,mean)^2))
        
        idxTableX[i,] <- c(ids_X,metric_X)
      }
  
      idxTableOmega <- idxTableOmega[order(idxTableOmega[,(self$cell_types+1)],decreasing=F),]
      idxTableX <- idxTableX[order(idxTableX[,(self$cell_types+1)],decreasing=F),]
      
      return(list(idsTableOmega = idxTableOmega[1:top,,drop=FALSE],
                  idsTableX = idxTableX[1:top,,drop=FALSE]))
    },
  
  ## Utils functions

    readInitValues = function(file) {
      initValues <- readRDS(file)
      ## Omega
      self$init_Omega <- initValues$init_Omega
    
      ## D
      self$init_D_w <- initValues$init_D_w
      self$init_D_h <- initValues$init_D_h
    
      ## X
      self$init_X <- initValues$init_X
    },

    readProjections = function(file) {
      initValues <- readRDS(file)
      ## R projection
      self$R <- initValues$R
      self$S <- initValues$S
      
      self$A <- matrix(apply(self$R,1,sum),ncol=1,nrow=self$cell_types)
      self$new_points <- self$V_row %*% t(self$R)
      self$B <- matrix(apply(self$S,1,sum),ncol=1,nrow=self$cell_types)
      self$new_samples_points <- t(self$S %*% self$V_column)
    },
    
    hinge = function(X) {
      sum(pmax(-X,0))
    },
    
    
    hinge_der_proportions = function(H,R,precision_=1e-10){
      m <- nrow(H)
      n <- ncol(H)
      der_R <- list()
      for (c in 1:m) {
        der_loc_R <- matrix(0,nrow=m,ncol=n)
        for (i in 1:m) {
          for (j in 1:n) {
            if (H[i,j] < -precision_) {
              der_loc_R[i,j] <- -R[c,j]
            }
          }
        }
        der_R[[c]] <- der_loc_R
      }
      res <- matrix(0,nrow=m,ncol=m)
      for (c in 1:m) {
        mtx <- der_R[[c]]
        res[,c] <- apply(mtx,1,sum)
      }
      res[,1] <- 0
      res
    },
    
    hinge_der_basis = function(W,S,precision_=1e-10){
      
      n <- ncol(W)
      res <- matrix(0,nrow=n,ncol=n)
      
      for (j in 1:n) {
        idx <- which(W[,j] < -precision_,arr.ind = T)
        res[,j] <- -apply(S[,idx,drop=F],1,sum)
      }
      res[1,] <- 0
      res
    },

    
  
    ## Optimization
  
    runGradientBlock = function(
      block_name = NULL,
      coef_der_X = self$coef_der_X,
      coef_der_Omega = self$coef_der_Omega,
      coef_hinge_H = self$coef_hinge_H,
      coef_hinge_W = self$coef_hinge_W,
      coef_pos_D_h = self$coef_pos_D_h,
      coef_pos_D_w = self$coef_pos_D_w,
      iterations = self$global_iterations,
      startWithInit = F,
      limit_X = 0,
      limit_Omega = 0,
      cosine_thresh = 0
    ) {

      R_limit_X <- 0
      R_limit_Omega <- 0
      
      if (is.null(self$X) | startWithInit) {
        self$blocks_statistics <- data.frame(matrix(0,nrow=0,ncol=13))
        self$errors_statistics <- NULL
        self$points_statistics_X <- NULL
        self$points_statistics_Omega <- NULL
        
        self$X <- self$init_X
        self$D_w <- self$init_D_w
        self$Omega <- self$init_Omega
        
        self$D_h <- self$init_D_h
        
        self$H_ <- self$X %*% self$R
        self$full_proportions <- diag(self$D_h[,1]) %*% self$H_
        self$init_count_neg_props <- sum(self$full_proportions < -1e-10)
        
        self$W_ <- t(self$S) %*% self$Omega 
        self$full_basis <- self$W_ %*% diag(self$D_w[,1])
        self$init_count_neg_basis <- sum(self$full_basis < -1e-10)
      }

      if (limit_X>0) {
        limit_num_X <- floor(nrow(self$V_row)*limit_X)
        R_limit_X <- norm(self$new_points[names(self$zero_distance_genes[limit_num_X]),-1],"2")
      }

      if (limit_Omega>0) {
        limit_num_Omega <- floor(ncol(self$V_row)*limit_Omega)
        R_limit_Omega <- norm(self$new_samples_points[names(self$zero_distance_samples[limit_num_Omega]),-1],"2")
      }
      
      step_errors_statistics <- matrix(0,nrow=iterations,ncol=10)
      step_points_statistics_X <- matrix(0,nrow=iterations,ncol=self$cell_types^2)
      step_points_statistics_Omega <- matrix(0,nrow=iterations,ncol=self$cell_types^2)
      
      if (is.null(block_name)) {
        block_name <- paste0("block_", nrow(self$blocks_statistics)+1)
      }
      
      if (is.null(self$errors_statistics)) {
        from_idx <- 1
      } else {
        from_idx <- nrow(self$errors_statistics) + 1
      }
      
      self$blocks_statistics <- rbind(self$blocks_statistics,
                                      c(block_name, from_idx, from_idx+iterations-1,
                                        coef_der_X, coef_der_Omega, coef_hinge_H,
                                        coef_hinge_W, coef_pos_D_h, coef_pos_D_w,
                                        iterations,limit_X,limit_Omega,cosine_thresh))
      
      colnames(self$blocks_statistics) <- c("block_name",
                                            "from", "to",
                                            "coef_der_X", "coef_der_Omega", 
                                            "coef_hinge_H", "coef_hinge_W", 
                                            "coef_pos_D_h", "coef_pos_D_w",
                                            "iterations", "limit_X", "limit_Omega",
                                            "cosine_thresh")
      
      res_ <- derivative_stage2(self$X, self$Omega, self$D_w,
                               self$V_row, self$R, self$S,
                               coef_der_X, coef_der_Omega,
                               coef_hinge_H, coef_hinge_W, coef_pos_D_h,
                               coef_pos_D_w, self$cell_types, self$N, self$M,
                               iterations, step_errors_statistics, 0,
                               step_points_statistics_X, step_points_statistics_Omega,
                               self$mean_radius_X, self$mean_radius_Omega,
                               R_limit_X, R_limit_Omega,cosine_thresh)
      
      self$X <- res_[[1]]
      self$Omega <- res_[[2]]
      self$D_w <- res_[[3]]
      self$D_h <- res_[[4]]
      self$errors_statistics <- rbind(self$errors_statistics,
                                      res_[[5]])
      self$points_statistics_X <- rbind(self$points_statistics_X,
                                        res_[[6]])
      self$points_statistics_Omega <- rbind(self$points_statistics_Omega,
                                            res_[[7]])
      
      colnames(self$errors_statistics) <- c("deconv_error","lamdba_error","beta_error",
                                            "D_h_error","D_w_error","total_error","orig_deconv_error",
                                            "neg_props_count","neg_basis_count","sum_d_w")
      self$H_ <- self$X %*% self$R
      self$full_proportions <- diag(self$D_h[,1]) %*% self$H_
      self$orig_full_proportions <- self$full_proportions
      self$count_neg_props <- sum(self$full_proportions < -1e-10)
      
      self$full_proportions[self$full_proportions < 0] <- 0
      self$full_proportions <- t(t(self$full_proportions) / (rowSums(t(self$full_proportions))+1e-10))
      self$full_proportions[is.nan(self$full_proportions)] <- 0
      
      
      self$W_ <- t(self$S) %*% self$Omega
      self$full_basis <- self$W_ %*% diag(self$D_w[,1])
      self$count_neg_basis <- sum(self$full_basis < -1e-10)
      
      self$orig_full_basis <- self$full_basis
      self$full_basis[self$full_basis < 0] <- 0
      self$full_basis <- self$full_basis / (rowSums(self$full_basis)+1e-10)
      self$full_basis[is.nan(self$full_basis)] <- 0
      
    },
  
  ## Plots
    
    plotErrors = function(variables = c("deconv_error","lamdba_error","beta_error",
    "D_h_error","D_w_error","total_error")) {
      
      toPlot <- data.frame(self$errors_statistics[,variables])
      toPlot$iteration <- 0:(nrow(self$errors_statistics)-1)
      toPlot <- melt(toPlot,id.vars="iteration",measure.vars = variables)
      plt <- ggplot(toPlot,aes(x=iteration,y=log10(value),color=variable)) +
        geom_point(size=0.2) +
        geom_line() + theme_minimal()
      plt
    },
  
  plotPoints2D = function(points="init",dims=3) {
    if (!points %in% c("init","current")) {
      stop("Allowed values for points are 'init', 'current'")
    }
    
    if (points == "init") {
      X <- self$init_X
      Omega <- self$init_Omega
      count_neg_props <- self$init_count_neg_props
      count_neg_basis <- self$init_count_neg_basis
    }
    if (points == "current") {
      X <- self$X
      Omega <- self$Omega
      count_neg_props <- self$count_neg_props
      count_neg_basis <- self$count_neg_basis
    }
    
    X <- X[,1:dims]
    Omega <- Omega[1:dims,]
    
    ## plot X
    toPlot <- as.data.frame(self$V_row %*% t(self$R))[,1:dims]
    colnames(toPlot) <- c("X","Y","Z")
    colnames(X) <- c("X","Y","Z")
    pltX <- ggplot(toPlot, aes(x=Y, y=Z)) +
      geom_point() + 
      geom_polygon(data=as.data.frame(X), fill=NA, color = "green") +
      theme_minimal()
    if (!is.null(count_neg_props)) {
      pltX <- pltX + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_props / (self$cell_types*self$N),4)*100,"%"), vjust=1, hjust=1)
    }
    
    ## plot Omega  
    toPlot <- as.data.frame(t(self$S %*% self$V_column))[,1:dims]
    colnames(toPlot) <- c("X","Y","Z")
    rownames(Omega) <- c("X","Y","Z")
    pltOmega <- ggplot(toPlot, aes(x=Y, y=Z)) +
      geom_point() + 
      geom_polygon(data=as.data.frame(t(Omega)), fill=NA, color = "green") +
      theme_minimal()
    
    if (!is.null(count_neg_basis)) {
      pltOmega <- pltOmega + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_basis / (self$cell_types*self$M),4)*100,"%"), vjust=1, hjust=1)
    }
    
    grid.arrange(pltX,pltOmega,nrow=1)
  },
  
  plotDistances = function() {
    if (is.null(self$distance_genes)) {
      stop("Run calculateDistances first")
    }
    
    toPlot_Genes <- data.frame(Distance=self$distance_genes)
    rownames(toPlot_Genes) <- names(self$distance_genes)
    toPlot_Genes$idx <- 1:nrow(toPlot_Genes)
    
    toPlot_Samples <- data.frame(Distance=self$distance_samples)
    rownames(toPlot_Samples) <- names(self$distance_samples)
    toPlot_Samples$idx <- 1:nrow(toPlot_Samples)
    
    pltGenes <- ggplot(toPlot_Genes,aes(x=idx,y=Distance)) +
      geom_point(size=0.1) +
      geom_line() + theme_minimal()
    
    pltSamples <- ggplot(toPlot_Samples,aes(x=idx,y=Distance)) +
      geom_point(size=0.1) +
      geom_line() + theme_minimal()
    
    grid.arrange(pltGenes,pltSamples)
    
  },
  
  plotMetricHistogram = function(logScale=F,
                                 breaks=100) {
    if (self$metric == "mean") {
      toPlot <- self$genes_mean
    } else if (self$metric == "sd") {
      toPlot <- self$genes_sd
    } else if (self$metric == "mad") {
      toPlot <- self$genes_mad
    } else {
      stop("Metric not find. Available metrics: mean, sd, mad")
    }
    
    if (logScale) {
      toPlot <- log(toPlot,10)
    }
    binwidth <- round((max(toPlot)-min(toPlot))/breaks)
    toPlot <- data.frame(toPlot)
    colnames(toPlot) <- "metric"
    p <- ggplot(toPlot, aes(x=metric)) + 
      geom_histogram(binwidth = binwidth)
    p
  },


    saveResults = function() {
## save proportions
      write.table(self$full_proportions,
                  file=paste0(self$path_,"/","markers_",self$analysis_name,"_proportions.tsv"),
                  sep="\t",col.names = NA, row.names = T, quote = F)
## save basis row normalized            
      colnames(self$full_basis) <- paste0("Cell_type_",1:self$cell_types)
      toSave <- self$full_basis
      toSave <- self$getFoldChange(toSave)
      #toSave <- self$getLimmaFoldChange(toSave)
      toSave <- rbind(c(rep(NA,self$cell_types),rep(NA,self$cell_types),round(apply(self$full_proportions,1,mean),4)),toSave)
      rownames(toSave) <- c("avg_proportions",rownames(self$filtered_dataset))
      write.table(toSave,file=paste0(self$path_,"/","markers_",self$analysis_name,"_basis_fc.tsv"),
                  sep="\t",col.names = NA, row.names = T, quote = F)
## save basis column normalized
      toSave <- t(t(self$full_basis) / rowSums(t(self$full_basis)))
      toSave <- self$getFoldChange(toSave)
      #toSave <- self$getLimmaFoldChange(toSave)
      toSave <- rbind(c(rep(NA,self$cell_types),rep(NA,self$cell_types),round(apply(self$full_proportions,1,mean),4)),toSave)
      rownames(toSave) <- c("avg_proportions",rownames(self$filtered_dataset))
      write.table(toSave,file=paste0(self$path_,"/","markers_",self$analysis_name,"_basis_fc_columns.tsv"),
                  sep="\t",col.names = NA, row.names = T, quote = F)  
    }
  )
)