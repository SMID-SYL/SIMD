code_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
data_dir = paste0('/',paste0(unlist(strsplit(code_dir,'\\/'))[2:6],collapse = "/"),'/data/ACPS/')
PSP_dir = paste0(data_dir,'PTE_PSP/')
drug_info_data_dir = paste0(data_dir,'Drug_information.txt')
#* cutoff = The significant P-value cutoff set by users
cutoff = 0.04

### Step1.calculate corrrelation using weight score
f_lst = list.files(PSP_dir)
dn_dat = read.csv(drug_info_data_dir,sep='\t',header=TRUE)  
for_tmp = c()
for(i in 1:length(f_lst)){
  dat = read.csv(paste0(PSP_dir,f_lst[i]),sep='\t',header=TRUE,row.names = 1)
  ws = apply(as.matrix(cbind(as.numeric(dat[,1]),as.numeric(dat[,3]))), 1, function(x) mean(x))
  cor_tst = cor.test(ws, as.numeric(dat[,2]))
  dn0 = strsplit(colnames(dat)[2],'\\.')[[1]]
  dn = paste(dn0[1],dn0[2],sep='-')
  dn2 = unique(dn_dat[which(dn_dat[,2] == dn),5])
  if(length(dn2) > 0){
    for(k in 1:length(dn2)){
      for_tmp = rbind(for_tmp, c(dn,dn2[k],cor_tst$estimate))
    }
  }
  print(i)
}
colnames(for_tmp) = c("DrugID","CID","cor_r")
write.table(for_tmp, paste0(data_dir,"step1_result_rPSPk.txt"),sep='\t',
            col.names = TRUE,row.names = FALSE,quote = FALSE)


### Step2.ACCP score calculation
acps_score <-function(x,dn_id,cutoff){
  library(DescTools)
  tmp0 = data.frame(drug_id = dn_id, cor_score = as.numeric(x) )
  x = tmp0$cor_score
  n = length(x)
  r = -(as.numeric(x))
  z = FisherZ(r)
  p = pnorm(z, mean = mean(z, na.rm = TRUE), sd = sd(z, na.rm = TRUE), lower.tail=FALSE)
  tmp = cbind(r,z,p)
  rownames(tmp) = tmp0[,1]
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}#scaling 0 to 1
  tmp = cbind(tmp[,1],range01(tmp[,2]),tmp[,3])
  tmp1 = tmp[which(tmp[,3] < cutoff),]
  tmp2 = aggregate(tmp1[,2], by=list(rownames(tmp1)), mean)
  tmp3 = tmp2[order(tmp2[,2],decreasing = TRUE),]
  pos = c(which(tmp3[,1]=="0"),which(tmp3[,1]=="-0"))
  if(length(pos) > 1){
    tmp3 = tmp3[-pos,] 
  }
  colnames(tmp3) = c("CID","accp_score")
  return(tmp3)
}


in_dat = read.csv(paste0(data_dir,"step1_result_rPSPk.txt"),sep='\t',header=TRUE)
id = unlist(lapply(in_dat[,2] , function(x) strsplit(sprintf("%f", x),'\\.')[[1]][1]))
ct = acps_score(in_dat[,3],id,i)
write.table(ct, paste0(data_dir,"step2_result_ACPS_",cutoff,".txt"),sep='\t',
            col.names = TRUE,row.names = FALSE,quote = FALSE)

### Step3.AUC calculation
library(ROCR)
accp_score_dat = read.csv( paste0(data_dir,"step2_result_ACPS_",cutoff,".txt"),sep='\t',header=TRUE )
GS_dat1 = read.csv(paste0(data_dir,'Prescribed_BRCA_TCGA.txt'),sep='\t',header = TRUE)
gs_tmp = rep(0,nrow(accp_score_dat))
pos = intersect(as.numeric(accp_score_dat[,1]), as.numeric(GS_dat1[,2]))
gs_tmp[which(as.numeric(accp_score_dat[,1]) %in% pos)] = 1
gs_tmp1 = cbind(accp_score_dat[,1],gs_tmp)
pred <- prediction(as.numeric(as.character(accp_score_dat[,2])),as.numeric(gs_tmp1[,2]))
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
perf <- performance(pred, measure = "auc")
print(perf@y.values[[1]])
  






