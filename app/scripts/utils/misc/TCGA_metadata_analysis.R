library(data.table)
library(tibble)
library(stats)
library(ggpubr)

clinical <- as.data.frame(fread("/Users/aladyeva.e/Downloads/clinical.project-TCGA-HNSC.2022-10-14/clinical.tsv",na.strings = "'--"))
samples <- as.data.frame(fread("/Users/aladyeva.e/Downloads/biospecimen.project-TCGA-HNSC.2022-10-14/sample.tsv",na.strings = "'--"))
proportions <- as.data.frame(fread("/Volumes/martyomov/Active/IndividualBackUps/aladyevae/deconvolution/data/results/HNSC_TPM_HK_CC_F_LIMITX_OMEGA_MAD_1_KNN_GENES_500_X/09_2022/ct16/20221012_171716/best/proportions.tsv")) %>% column_to_rownames("V1")
rownames(proportions) <- paste0("Cell_type_",rownames(proportions))

proportions_t <- as.data.frame(round(t(proportions),4))
proportions_t <- proportions_t[,apply(proportions_t,2,mean)>0.01]
proportions_t <- proportions_t %>% rownames_to_column("full_id")


clinical_mod <- dcast(clinical,case_submitter_id ~ treatment_type,value.var="treatment_or_therapy")
clinical <- clinical[,-c(which(grepl("^treatment_",colnames(clinical))))]
clinical <- clinical[duplicated(clinical),]
clinical <- merge(clinical,clinical_mod,on="case_submitter_id")
clinical <- clinical[,colSums(is.na(clinical))!=nrow(clinical)]

data_ <- data.frame(full_id = colnames(proportions))
data_$case_submitter_id <- sapply(data_$full_id,function(x){paste0(strsplit(x,"-")[[1]][1:3],collapse="-")})
data_$sample_submitter_id <- sapply(data_$full_id,function(x){paste0(strsplit(x,"-")[[1]][1:4],collapse="-")})

data_ <- merge(data_,clinical,on="case_submitter_id")
data_ <- merge(data_,proportions_t,on="full_id")
data_ <- merge(data_,samples[,c("sample_submitter_id","sample_type")],on="sample_submitter_id")

cols_ <- c("vital_status","gender","race","ajcc_clinical_m","ajcc_clinical_n","ajcc_clinical_stage",
           "ajcc_clinical_t","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stage","ajcc_pathologic_t",
           "ajcc_staging_system_edition","primary_diagnosis","prior_malignancy","prior_treatment","synchronous_malignancy",
           "tissue_or_organ_of_origin","Pharmaceutical Therapy, NOS","Radiation Therapy, NOS","sample_type")
cols_ <- c("vital_status","sample_type")
for (col_ in cols_) {
  print(col_)
  #short_data_ <- data_[,c(colnames(data_)[grepl("^Cell_",colnames(data_))],col_)]
  short_data_ <- data_[,c(col_,"Cell_type_2","Cell_type_6","Cell_type_7","Cell_type_9","Cell_type_15","Cell_type_5")]
  cor_data <- model.matrix(~0+., data=short_data_) %>% cor(use="complete.obs") 
  toPlot <- melt(cor_data[grepl("^Cell_",colnames(cor_data)),grepl(paste0("^",col_),rownames(cor_data))])
  toPlot$Var2 <- gsub(paste0("^",col_),"",toPlot$Var2)
  plt <- ggplot(toPlot) + geom_tile(aes(x=Var2,y=Var1,fill=value),
                             color = "white") +
    geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 3) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1)) + theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                     size = 12, hjust = 1)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank()) +
    coord_fixed()  
  print(plt)
}

cols_ <- c("vital_status","sample_type")
col_ <- "vital_status"
for (col_ in cols_) {
toPlot <- data_[,c(col_,"Cell_type_2","Cell_type_6","Cell_type_7","Cell_type_9","Cell_type_15","Cell_type_5")]
toPlot <- melt(toPlot)
vars_ <- length(toPlot[,col_])
toPlot[,col_] <- as.factor(toPlot[,col_])
toPlot$variable <- as.factor(toPlot$variable)
if (vars_>2) {
  p <- ggviolin(data = toPlot, x=col_, y="value", color = "variable", 
                fill="variable", palette = "npg", 
                add = "mean_sd", add.params = list(color = "black", size=0.2))
  p <- facet(p + theme_bw(), facet.by = "variable",
             nrow=3, strip.position = "right") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    stat_compare_means(method = "anova",label = "p.signif")
} else {
  p <- ggviolin(data = toPlot, x=col_, y="value", color = "variable", 
                fill="variable", palette = "npg", 
                add = "mean_sd", add.params = list(color = "black", size=0.2))
  p <- facet(p + theme_bw(), facet.by = "variable",
             nrow=3, strip.position = "right") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    stat_compare_means(label = "p.signif", method = "wilcox.test")
}
print(p)
}






