library(dplyr)
library(readr)
library(tools)
library(data.table)

source("/app/scripts/OutputPlots.R")

df <- snakemake@input$stats %>% 
  lapply(read_tsv) %>% 
  bind_rows %>% mutate(analysis_name = gsub("_stats$","",file_path_sans_ext(basename(snakemake@input$stats))))

df$init_number <- sapply(df$analysis_name,function(x){tt <- strsplit(x,"_")[[1]]; tt[length(tt)]})
rest_columns <- c("init_number","analysis_name",colnames(df)[-c(ncol(df),ncol(df)-1)])
df <- df[,rest_columns, with=FALSE]
write.table(df,file=snakemake@output[[1]],
                  sep="\t",col.names = T, row.names = F, quote = F)  

df <- df %>%
  slice(which.min(total_error))

idx <- which(file_path_sans_ext(basename(snakemake@input$stats)) == paste0(df$analysis_name,"_stats"))
cat(df$analysis_name,file=snakemake@output[[3]],sep="\n")

metadata_ <- loadRData(file.path("results",snakemake@params[["run_id"]],"meta",paste0(df$analysis_name,".meta")))

png(snakemake@output[[4]])
plotNegBasisChange(metadata_)
dev.off()

png(snakemake@output[[5]])
plotNegProportionsChange(metadata_)
dev.off()

png(snakemake@output[[6]])
plotSumToOneChange(metadata_)
dev.off()

png(snakemake@output[[7]])
plotPoints2D(metadata_,"init")
dev.off()

png(snakemake@output[[8]])
plotPoints2D(metadata_,"current")
dev.off()

file.copy(snakemake@input[[idx]], 
snakemake@output[[2]])

file.copy(file.path("results",snakemake@params[["run_id"]],"meta",paste0(df$analysis_name,".meta")),
          snakemake@output[[9]])

file.copy(file.path("results",snakemake@params[["run_id"]],"props",paste0(df$analysis_name,"_proportions.tsv")),
          snakemake@output[[10]])

file.copy(file.path("results",snakemake@params[["run_id"]],"basis_row",paste0(df$analysis_name,"_basis_fc.tsv")),
          snakemake@output[[11]])

file.copy(file.path("results",snakemake@params[["run_id"]],"basis_col",paste0(df$analysis_name,"_basis_fc_clmn.tsv")),
          snakemake@output[[12]])

file.copy(file.path("results",snakemake@params[["run_id"]],"points",paste0(df$analysis_name,"_points.rds")),
snakemake@output[[14]])

all_basis  <- data.frame(matrix(0,ncol=0,nrow=nrow(metadata_$filtered_dataset)))
for (f in snakemake@input[['basis_col']]){
  tmp_ <- as.data.frame(fread(f))
  tmp_ <- tmp_[-1,]
  tmp_ <- tmp_[,colnames(tmp_)[grepl("^Cell",colnames(tmp_))]]
  analysis_name <- gsub("_basis_fc_clmn","",file_path_sans_ext(basename(f)))
  colnames(tmp_) <- paste0(analysis_name,"_",colnames(tmp_))
  all_basis <- cbind(all_basis,tmp_)
}

all_X <- data.frame(matrix(0,nrow=0,ncol=metadata_$cell_types))
all_Omega <- data.frame(matrix(0,nrow=0,ncol=metadata_$cell_types))
  
for (f in snakemake@input[['points']]) {
  points_ <- readRDS(f)
  analysis_name <- gsub("_points","",file_path_sans_ext(basename(f)))
  tmp_X <- as.data.frame(points_$final_X)
  tmp_Omega <- as.data.frame(t(points_$final_Omega))
  rownames(tmp_X) <- paste0(analysis_name,"_Cell_type_",rownames(tmp_X))
  rownames(tmp_Omega) <- paste0(analysis_name,"_Cell_type_",rownames(tmp_Omega))
  all_X <- rbind(all_X,tmp_X)
  all_Omega <- rbind(all_Omega,tmp_Omega)
}

png(snakemake@output[['UMAP']])
plotUMAP(t(all_basis),df$analysis_name)
dev.off()

png(snakemake@output[['UMAP_projX']])
plotUMAP(all_X,df$analysis_name)
dev.off()

png(snakemake@output[['UMAP_projOmega']])
plotUMAP(all_Omega,df$analysis_name)
dev.off()

png(snakemake@output[[13]])
plotAbundance(metadata_)
dev.off()

png(snakemake@output[[15]])
plotDistribution(metadata_$orig_full_proportions)
dev.off()

png(snakemake@output[[16]])
plotDistribution(metadata_$orig_full_basis)
dev.off()

basis <- fread(file.path("results",snakemake@params[["run_id"]],"basis_col",paste0(df$analysis_name,"_basis_fc_clmn.tsv")))
basis <- basis[-1,]

png(snakemake@output[[17]])
plotBasisDistance(metadata_,basis)
dev.off()

proportions <- fread(file.path("results",snakemake@params[["run_id"]],"props",paste0(df$analysis_name,"_proportions.tsv")))

png(snakemake@output[[18]])
plotProportionsDistance(metadata_,proportions)
dev.off()