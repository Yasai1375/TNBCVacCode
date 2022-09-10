library(GenomicDataCommons)
library(TCGAbiolinks)
library(remotes)
library(apeglm)
library(ashr)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)



# Receiving TNBC clinical data --------------------------------------------

query_clinical<- GDCquery(project = "TCGA-BRCA", 
                          data.category = "Clinical",
                          data.type = "Clinical Supplement", 
                          data.format = "BCR Biotab",
                          file.type="patient")
GDCdownload(query_clinical)
clinical.BCRtab.all <- GDCprepare(query_clinical)

clinical<-clinical.BCRtab.all$clinical_patient_brca



# ER Negative and PR Negative and her2_status_by_ihc Negative 

er_pr_her2<-clinical %>%
  filter(er_status_by_ihc=="Negative" & pr_status_by_ihc=="Negative" & her2_fish_status!="Positive" &
           her2_status_by_ihc=="Negative" & her2_ihc_score %in% c ("0", "1+", "[Not Available]"))



# ER Negative and PR Negative and her2_status_by_ihc Equivocal or Indeterminate but her2_fish_status Negative

retrieve_er_pr_her2<-clinical %>%
  filter(er_status_by_ihc=="Negative" & pr_status_by_ihc=="Negative" & her2_fish_status=="Negative" &
           her2_status_by_ihc %in% c("Negative","Equivocal","Indeterminate","[Not Evaluated]"))



# barcode_list_of_TNBC_cases ----------------------------------------------

barcode<-unique(c(er_pr_her2$bcr_patient_barcode , retrieve_er_pr_her2$bcr_patient_barcode))



# obtaining HTSeq - Counts from GDC portal --------------------------------

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  barcode = barcode)
GDCdownload(query)
data <- GDCprepare(query)



# Gene Expression Analysis ------------------------------------------------

BRCAMatrix <- assay(data,"HTSeq - Counts")
ddsSE <- DESeqDataSet(data, design = ~ shortLetterCode)
keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)
res<-results(ddsSE, name = "shortLetterCode_TP_vs_NT" , alpha = 0.05)
df_TNBC_DEG <- as.data.frame(res)
summary(res)
save(df_TNBC_DEG, file = "df_TNBC_DEG.Rda")



# shrink with apeglm method -----------------------------------------------

res_shrink_apeglm <- lfcShrink(dds = ddsSE , coef = "shortLetterCode_TP_vs_NT" , type = "apeglm" )
df_TNBC_DEG_shrink_apeglm <- as.data.frame(res_shrink_apeglm)



# list of ensembl_id of selected cancer testis antigens -------------------

selected_id_rep<-c("ENSG00000130377","ENSG00000111254","ENSG00000137948","ENSG00000141371",
                   "ENSG00000137225","ENSG00000205111","ENSG00000176566","ENSG00000152670",
                   "ENSG00000268606","ENSG00000142025","ENSG00000135436","ENSG00000189132",
                   "ENSG00000221867","ENSG00000143194","ENSG00000147381","ENSG00000124260",
                   "ENSG00000099399","ENSG00000189023","ENSG00000156269","ENSG00000143552",
                   "ENSG00000177947","ENSG00000163114","ENSG00000175646","ENSG00000170748",
                   "ENSG00000165496","ENSG00000181433","ENSG00000112053","ENSG00000039600",
                   "ENSG00000123569","ENSG00000139351","ENSG00000102387","ENSG00000180113",
                   "ENSG00000006047","ENSG00000227234","ENSG00000184033","ENSG00000046774",
                   "ENSG00000155495","ENSG00000185247","ENSG00000198681","ENSG00000137090")



# CTA selection based on LFC and padj -------------------------------------

final_selected_rep<-df_TNBC_DEG_shrink_apeglm %>%
  filter(row.names(df_TNBC_DEG_shrink_apeglm) %in% selected_id_rep) %>%
  filter(log2FoldChange>1 & padj<=0.05) %>%
  arrange(desc(log2FoldChange))



# selecting expressed samples ---------------------------------------------

# 1. obtaining HTSeq - FPKM from GDC portal -------------------------------

query_fpkm <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "HTSeq - FPKM",
                       barcode = barcode)
GDCdownload(query_fpkm)
data_fpkm <- GDCprepare(query_fpkm)
BRCAMatrix_fpkm <- assay(data_fpkm,"HTSeq - FPKM")

# 2. filter TP data --------------------------------------------------------

data_fpkm_TP<-data_fpkm[,which(data_fpkm$shortLetterCode=="TP")]
BRCAMatrix_fpkm_TP <- assay(data_fpkm_TP,"HTSeq - FPKM")

data_TP<-data[,which(data$shortLetterCode=="TP")]
BRCAMatrix_TP <- assay(data_TP,"HTSeq - Counts")

# 3. the percentage of patients who express each cancer testis gene -------

expressed_percent<-function(id,treshold_count=10,treshold_fpkm=0.3){
  percentage=sum(BRCAMatrix_TP[id,]>=treshold_count & BRCAMatrix_fpkm_TP[id,]>=treshold_fpkm)/154*100 
  return(percentage)
}



# selecting CTA with at least 10% expression in patients population -------

df_percent_rep<-data.frame(id=row.names(final_selected_rep),percentage=unlist(lapply(row.names(final_selected_rep), expressed_percent)))
df_percent_rep<-df_percent_rep %>%
  filter(percentage>=10)



# list of expressed patients for each CTA  --------------------------------

expressed_patient<-function(id,treshold_count=10,treshold_fpkm=0.3){
  names<-colnames(BRCAMatrix_TP[,BRCAMatrix_TP[id,]>=treshold_count & BRCAMatrix_fpkm_TP[id,]>=treshold_fpkm])
  return(data.frame(id,names))
}



# calculating the total coverage of vaccine in TNBC population ------------

expressed_names<-c(expressed_patient("ENSG00000221867")$names,expressed_patient("ENSG00000046774")$names,
                   expressed_patient("ENSG00000137090")$names,expressed_patient("ENSG00000147381")$names,
                   expressed_patient("ENSG00000102387")$names,expressed_patient("ENSG00000143194")$names)

unique(expressed_names)
