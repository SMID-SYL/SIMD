rm(list=ls())
library(Seurat)
library(miQC)
library(SeuratWrappers)
library(harmony)
library(garnett)
library(org.Hs.eg.db)
library(dplyr)
library(fpc)
library(cluster )
library(naturalsort)
library(infercnv)

code_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir = paste0('/',paste0(unlist(strsplit(code_dir,'\\/'))[2:6],collapse = "/"),'/data/PPNE/')
rds_dir = paste0(data_dir,'1_rds/' )
func_dir = paste0(data_dir,'2_infer_CNV_results/2_NE_function_file/')

####################################################################################################
#EMTAB8107 data analysis
####################################################################################################
#Step1. QC process
in.obj = readRDS( paste0(rds_dir,"0_merge_dat.rds"))
data.fit = in.obj[["percent.mt"]] = PercentageFeatureSet(object = in.obj, pattern = "^MT-")
data.fit = RunMiQC(in.obj, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75, model.slot = "flexmix_model")
in.obj = subset(data.fit, miQC.keep == "keep")

#Step2 Data processing
in.obj = NormalizeData(object = in.obj, normalization.method = "LogNormalize")
in.obj = FindVariableFeatures(object = in.obj, selection.method = "vst")
in.obj= ScaleData(object = in.obj, features = rownames(in.obj))
in.obj = RunPCA(in.obj, verbose = F, npcs = 20)
in.obj = JackStraw(object = in.obj, num.replicate = 100, dims = 50)
in.obj = ScoreJackStraw(object = in.obj, dims = 1:20)

#Step3 batch correction + Clustring + diemsion reduction
rds_dir = '/home/starjjbang/2022/Single_cell/1_data_analysis/3_rda/'
harmony.obj <- RunHarmony(object = in.obj, group.by.vars = "orig.ident")
harmony.obj <- FindNeighbors(object = harmony.obj, dims = 1:20, force.recalc = T, reduction = "harmony")
harmony.obj <- FindClusters(object = harmony.obj, resolution = 0.6)
harmony.obj <- RunUMAP(object = harmony.obj, reduction = "harmony", dims = 1:20)
harmony.obj <- RunTSNE(object = harmony.obj, reduction = "harmony", dims = 1:20)

#Step4 cell annotation using marker genes with garnett
temp_raw_data <- GetAssayData(object = harmony.obj, slot = "counts")
pd <- new("AnnotatedDataFrame", data = harmony.obj@meta.data)
fData <- data.frame(gene_short_name = row.names(temp_raw_data), row.names = row.names(temp_raw_data))
fd <- new("AnnotatedDataFrame", data = fData)
lowerDetectionLimit <- 0
if(all(temp_raw_data == floor(temp_raw_data))) {
  expressionFamily <- negbinomial.size()
} else if(any(data < 0)){
  expressionFamily <- uninormal()
} else {
  expressionFamily <- tobit()
}
monocle_cds <- newCellDataSet(temp_raw_data,
                              phenoData = pd, 
                              featureData = fd,
                              lowerDetectionLimit=lowerDetectionLimit,
                              expressionFamily=expressionFamily)
monocle_cds <- estimateSizeFactors(monocle_cds)
marker_file_path <- paste0(data_dir,'cell_type_markers.txt')
marker_check <- check_markers(monocle_cds, 
                              marker_file_path,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL",
                              classifier_gene_id_type= "SYMBOL")
temp_classifier <- train_cell_classifier(cds = monocle_cds,
                                         marker_file = marker_file_path,
                                         db= org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         marker_file_gene_id_type = "SYMBOL",
                                         classifier_gene_id_type = "SYMBOL",
                                         cores = 25,
                                         num_unknown=5000,
                                         max_training_samples=5000)
monocle_cds <- classify_cells(monocle_cds, temp_classifier,
                              db = org.Hs.eg.db,
                              cluster_extend = TRUE,
                              cds_gene_id_type = "SYMBOL")
head(pData(monocle_cds))
harmony.obj@meta.data$gcell_type <- pData(monocle_cds)$cell_type
saveRDS(harmony.obj, paste0(rds_dir,"2_Processed_EMTAB8107.rds"))

############################################################################
#Extract NE cells
############################################################################ 
### inferCNV : input_generation  , Epithelial only
#harmony.obj = readRDS(paste0(rds_dir,"2_Processed_EMTAB8107.rds"))
cnv_re_dir = paste0(data_dir,'2_infer_CNV_results/')

###  Epithelial only ########################################
#input data make
ext.epi <- subset(harmony.obj, subset = gcell_type == c("Epithelial"))
ext.epi1 <- subset(harmony.obj, subset = gcell_type == c("T cells","B cells"))
expr <- as.data.frame(GetAssayData(object = ext.epi, assay.type = "RNA", slot = "data"))
expr1 <- as.data.frame(GetAssayData(object = ext.epi1, assay.type = "RNA", slot = "data"))
m_dat2 = merge(expr,expr1,by="row.names")
m_dat3 = m_dat2[,2:ncol(m_dat2)]
rownames(m_dat3) = m_dat2[,1]
write.table(x = m_dat3,file = paste0(cnv_re_dir,'1_input/ct_exp.csv'), quote = F, sep = "\t")
an_f = data.frame(id = rownames(ext.epi@meta.data),cluster = ext.epi@meta.data$gcell_type)
an_f1 = data.frame(id = rownames(ext.epi1@meta.data),cluster = ext.epi1@meta.data$gcell_type)
an_f2 = rbind(an_f,an_f1)
write.table(an_f2, file = paste0(cnv_re_dir,'1_input/ct_annotation.csv'), quote = F, sep = "\t", col.names = F, row.names = F)

### object make ################
#1. only epi
cnv_re_dir = paste0(data_dir,'2_infer_CNV_results/')
gene.order.file <-paste0(cnv_re_dir,"1_input/hg19.inferCNV.gtf")
infercnv.obj <- CreateInfercnvObject(
  raw_counts_matrix = paste0(cnv_re_dir,'1_input/ct_exp.csv'),
  annotations_file = paste0(cnv_re_dir,'1_input/ct_annotation.csv'),
  gene_order_file = gene.order.file,
  ref_group_names = c("T cells","B cells"),delim = "\t")

### running ################
dir.create(file.path(cnv_re_dir, 'infer_CNV_RUN_result_1','/'))
infercnv_obj_default <- infercnv::run(infercnv.obj,cutoff=0.1, 
                                      out_dir = paste0(cnv_re_dir, 'infer_CNV_RUN_result_1','/'),
                                      window_length=101,
                                      max_centered_threshold=3,
                                      cluster_by_groups=F,
                                      plot_steps=FALSE,
                                      denoise=TRUE,
                                      sd_amplifier=1.3,
                                      analysis_mode = "samples",
                                      HMM=FALSE,
                                      png_res=300,
                                      num_threads=80,
                                      write_expr_matrix=TRUE)
### ture epithelial cell filtering ################
include_group_annotation <- FALSE
adjust_normal_thresholds <- TRUE
cancer_x_threshold_sd_multiplier <- 2
cancer_y_threshold_sd_multiplier <- 1.5
normal_x_threshold_sd_multiplier <- 1
normal_y_threshold_sd_multiplier <- 1.25

### 0. Define functions ###
gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir,  "fetch_chromosome_boundaries.R"))
create_array_CNV_annotation <- dget(paste0(func_dir, "create_array_CNV_annotation.R"))

### 1. Load InferCNV output and create heatmap and metadata dfs ###
cnv_re_dir = paste0(data_dir,'2_infer_CNV_results/')
infercnv_output <- t(read.csv(paste0(cnv_re_dir,'infer_CNV_RUN_result/infercnv.observations.txt'),sep=' '))
harmony.obj = readRDS(paste0(rds_dir,"2_Processed_EMTAB8107.rds"))
metadata <- prepare_infercnv_metadata(harmony.obj, subset_data = F, as.data.frame(t(infercnv_output)), for_infercnv=F)
epithelial_metadata <- metadata$metadata
epithelial_ids <- epithelial_metadata$cell_ids[grep("pithelial", epithelial_metadata$cell_type)]
epithelial_heatmap <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]

### 2. Add QC metadata ###
# create epithelial_metadata df and only include epithelial cells in epithelial_heatmap:
epithelial_metadata <- epithelial_metadata[rownames(epithelial_heatmap),]
# order epithelial metadata cell type cluster levels:
epithelial_metadata$cell_type <- factor(
  epithelial_metadata$cell_type,
  levels = naturalsort(unique(epithelial_metadata$cell_type))
)
QC <- data.frame(
  row.names = names(Idents(harmony.obj)),
  nUMI = harmony.obj@meta.data$nCount_RNA,
  nGene = harmony.obj@meta.data$nFeature_RNA
)
QC <- QC[rownames(epithelial_metadata),]
epithelial_metadata <- cbind(epithelial_metadata, QC)

### 3. Add CNA and correlation value metadata ###
# determine CNA values and add to epithelial_metadata:
#print("Determining CNA values and adding to epithelial metadata df...")
# scale infercnv values to -1:1, square values and take the mean:
scaled_df <- as.data.frame(rescale(as.matrix(epithelial_heatmap), c(-1,1)))

#the mean of the squares of these values was used to define a genomic instability score for each cell.
CNA_values <- apply(scaled_df, 1, function(y) {
  #y[is.na(y)] <- 0
  #scaled_y <- rescale(y, c(-1, 1))
  return(mean(y^2))
})
CNA_value_df <- data.frame(row.names = names(CNA_values),CNA_value = CNA_values)
epithelial_metadata <- cbind(epithelial_metadata, CNA_value_df)
head(epithelial_metadata)

# determine top 5% cancer cells:
CNA_order <- order(CNA_value_df$CNA_value, decreasing=T)
ordered_CNA_values  <- data.frame(row.names = rownames(CNA_value_df)[CNA_order],CNA_value = CNA_value_df[CNA_order,])
top_cancer <- head(ordered_CNA_values, nrow(ordered_CNA_values)*0.05)
# find average genome-wide CNV predictions across genome:
top_cancer_CNV_average <- apply(epithelial_heatmap[rownames(top_cancer),], 2, mean)

# find correlations of each cell's CNVs with top_GIN_CNV_average:
for(i in 1:nrow(epithelial_heatmap)){
  x = epithelial_heatmap[i,]
  if (length(unique(as.numeric(x))) == 1) {
    cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
                             cor.p.value="no_CNVs_recorded")
  } else {
    cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "kendall")
    cor_result <- data.frame(cor$estimate, cor$p.value)
  }
  rownames(cor_result)=rownames(epithelial_heatmap)[i]
  write.table(cor_result[1,],paste0(cnv_re_dir,'cor.txt'),col.names = FALSE,quote=FALSE,sep='\t', append=TRUE)
  print(i)
}

### Call normals and add to epithelial_metadata ###
correlation_df <-unique(read.csv(paste0(cnv_re_dir,'cor.txt'),sep='\t'))
colnames(correlation_df) = c('cell_ids','cor.estimate','cor.p.value')
epithelial_metadata <- merge(epithelial_metadata,correlation_df,by="cell_ids")
head(epithelial_metadata)

#epithelial_metadata = epithelial_metadata[,-5]
quad_df <- data.frame(
  row.names = rownames(epithelial_metadata),
  CNA_value = epithelial_metadata$CNA_value, 
  cor.estimate = epithelial_metadata$cor.estimate
)
scaled_quad_df <- scale(quad_df) %>% as.data.frame()

# run silhouette cluster analysis to determine clusters and thresholds:
pamk_result <- pamk(scaled_quad_df, krange=1:4)
pamk_result$nc
silhouette_result <- pam(scaled_quad_df, pamk_result$nc)
sil_values <- as.data.frame(silhouette_result$silinfo$widths)
sil_result <- data.frame(row.names=names(silhouette_result$clustering),
                         cluster=silhouette_result$clustering,
                         sil_width=sil_values$sil_width)

# add sil_result to epithelial_metadata:
epithelial_metadata <- cbind(epithelial_metadata, sil_result)

# determine normal and cancer clusters by determining the max CNA values and
# correlation with top 5% cancer:
cluster_split <- split(epithelial_metadata, epithelial_metadata$cluster)
names(cluster_split) <- paste0("cluster_", names(cluster_split))
# determine order of clusters by adding mean CNA and correlation values:
cluster_means <- sort(
  unlist(
    lapply(cluster_split, function(x) {
      return(mean(x$CNA_value) + mean(x$cor.estimate))
    })
  )
)
# determine second cluster from axes as cancer cluster closest to axes:
first_cancer_cluster <- names(cluster_means[2])
first_cancer_df <- eval(parse(text=paste0("cluster_split$", first_cancer_cluster)))

# make intercepts 1 std dev from mean towards axes:
# define x-axis as 2 std devs left of mean:
cancer_x_threshold_sd_multiplier <- 2.5
cancer_y_threshold_sd_multiplier <- 1.5
normal_x_threshold_sd_multiplier <- 1
normal_y_threshold_sd_multiplier <- 1.25
CNA_mean <- mean(first_cancer_df$CNA_value)
CNA_std_dev <- sd(first_cancer_df$CNA_value)
x_int <- round(CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev), 3)
# define x-axis as 2 std devs left of mean:
cor_mean <- mean(first_cancer_df$cor.estimate)
cor_std_dev <- sd(first_cancer_df$cor.estimate)
y_int <- round(cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev), 3)
epithelial_metadata$normal_cell_call <- "cancer"
epithelial_metadata$normal_cell_call[epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int] <- "normal"
epithelial_metadata$normal_cell_call[epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate > y_int] <- "unassigned"
table(epithelial_metadata$normal_cell_call)
write.table(epithelial_metadata,paste0(cnv_re_dir,'epithelial_metadata.txt'),col.names = TRUE,quote=FALSE,sep='\t')

ct_lab = harmony.obj@meta.data$gcell_type
for(i in 1:nrow(epithelial_metadata)){
  ct = harmony.obj@meta.data[which(rownames(harmony.obj@meta.data) == epithelial_metadata$cell_ids[i]),10] 
  if(ct=="Epithelial"){
    ct_lab[which(rownames(harmony.obj@meta.data) == epithelial_metadata$cell_ids[i])] <- paste0(ct,'-',epithelial_metadata[i,10])
  }
}
harmony.obj@meta.data <-cbind(harmony.obj@meta.data,ct_lab)
colnames(harmony.obj@meta.data)[which(names(harmony.obj@meta.data) == "ct_lab")] <- "cell_type_cancer"
saveRDS(epithelial_metadata,paste0(rds_dir,"3_NE_cell.rds"))


############################################################################
#Calcaulated correlation between NE cells and PSP
############################################################################ 
library(biomaRt)
library(tidyverse)
library(dbplyr)
library(Kendall)

harmony.obj = readRDS(paste0(rds_dir,"2_Processed_EMTAB8107.rds"))
epithelial_metadata = readRDS(paste0(rds_dir,"3_NE_cell.rds"))
#Epithelial-cancer --> NE cells
#Epithelial-normal --> NNE cells
ct_ty =  "Epithelial-cancer"

Idents(object = harmony.obj) <- harmony.obj@meta.data$cell_type_cancer
deg.markers1 <- FindMarkers(harmony.obj, ident.1 = ct_ty, ident.2 = NULL, only.pos = TRUE,test.use="bimod", verbose = FALSE)
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
a = getBM(attributes=c("hgnc_symbol","ensembl_gene_id","entrezgene_id"), filters=c("hgnc_symbol"), 
          mart=ensembl, rownames(deg.markers1))
en_id = sapply(rownames(deg.markers1), function(x) {ifelse(length(a[which(a$hgnc_symbol == x),2]  ) > 0,  a[which(a$hgnc_symbol == x),2], NA)})
deg.markers2  = cbind(deg.markers1, en_id)
head(deg.markers2)

b_dat1 <- aggregate(avg_log2FC ~ en_id, data = deg.markers2, FUN = mean)
b_dat = data.frame(avg_log2FC = b_dat1[,2])
rownames(b_dat) = b_dat1[,1]
database_dir = paste('/DEE',sep='')
d_db = list.files(database_dir)
re_dir = paste0(data_dir,'3_COR/')
for_tmp = c()
for (j in 2:length(d_db)){
  df = as.matrix(read.table(paste0('/home/starjjbang/2022/2021_breast_cancer_database/',d_db[j])),sep='\t')
  fn = gsub('.txt','',paste(strsplit(d_db[j],'_')[[1]][3:6],collapse  ='_'))
  m_dat = merge(b_dat,df,by="row.names")
  t_val = Kendall(m_dat$avg_log2FC,m_dat$V1)
  for_tmp = rbind(for_tmp,c(fn, t_val$tau,  t_val$sl[1])  )
}
colnames(for_tmp) = c("drug","tau","P-value")
#Epithelial-cancer --> NE cells
#Epithelial-normal --> NNE cells
write.table(for_tmp, paste0(re_dir,'NEcell.txt'),sep='\t',quote=FALSE,col.names = TRUE,row.names=FALSE)


