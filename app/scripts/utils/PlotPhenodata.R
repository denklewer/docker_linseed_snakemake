library(data.table)
library(reshape2)
library(tibble)
library(ggpubr)

phenodata <- as.data.frame(fread("/Users/aladyeva.e/Dropbox (ArtyomovLab)/GTEx_lung/GTEx_sample_attributes_and_subject_lung.tsv"))
attribute <- "AGE"
sample_column <-"SAMPID"
PATH <- "/Volumes/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/results/GTEX_LUNGS_TPM_HK_CC_F_LIMITX_OMEGA_MAD_1_KNN_Rand_new/09_2022"
DT <- "20221007_220240"
cell_types <- 14

plotPhenodata <- function(PATH,DT,cell_types,phenodata,attribute,
                          sample_column="SAMPID",filter_=T,filter_limit=0.01){
  DATA_PATH <- file.path(PATH,paste0("ct",cell_types),DT,"best","proportions.tsv")
  props_ <- as.matrix(fread(DATA_PATH) %>% column_to_rownames("V1"))
  if (filter_) {
    avg_props <- apply(props_,1,mean)
    props_ <- props_[avg_props>filter_limit,]
  }
  toPlot <- melt(props_)
  colnames(toPlot) <- c("Cell_type","Sample","Proportions")
  phenodata[,sample_column] <- gsub("\\.","-",phenodata[,sample_column])
  toPlot <- merge(x=toPlot,y=phenodata,by.x=c("Sample"),by.y = c(sample_column),all.x=TRUE)
  toPlot <- toPlot[,c("Sample","Cell_type","Proportions",attribute)]
  toPlot[,attribute] <- as.factor(toPlot[,attribute])
  toPlot$Cell_type <- as.factor(toPlot$Cell_type)
  
  p <- ggviolin(data = toPlot, x=attribute, y="Proportions", color = "Cell_type", 
                fill="Cell_type", palette = "npg", 
                add = "mean_sd", add.params = list(color = "black", size=0.2))
  p <- facet(p + theme_bw(), facet.by = "Cell_type",
             nrow=3, strip.position = "right") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    stat_compare_means(method = "anova",label = "p.signif")
  p
}

plotPhenodata(PATH,DT,cell_types, phenodata, attribute)
