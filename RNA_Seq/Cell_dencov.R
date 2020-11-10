###################################
#
#
#
#
#
###################################
library(tximport)
library(EnsDb.Hsapiens.v86)
txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
gdf=read_tsv("~/dn/Projects/LTF-Shashwat/geneid")
tx2gene <-  as.data.frame(txdf[,c("tx_id","gene_id")])


setwd("~/Documents/Aldrin/CPDR_RNA-seq/kallisto_quant/")
files=Sys.glob("*/*.h5")
files1=data.frame(file=sapply(strsplit(files,"_"),'[',1))

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,ignoreTxVersion = T)

tpm_dt=txi$abundance
tpm_xcell=data.frame(tpm_dt) %>% rownames_to_column(var = "gene_id") %>%
  left_join(.,gdf,by=c("gene_id")) %>%
  select(gene_name,X1:X18) %>%
  write_tsv(.,"tpm_xcell.tsv")


###ID whether 
library(readxl)
covar=read_xlsx("~/Documents/G/Info_Mgmt/Bioinformatics-CPDR/Datasets/01_7AA_7CA_CaP_RNASeq/02 201210 ExpressionAnalytics RNASeq/EA RNA-Seq QC Data 101812/NextGen Seq RNA Specimens Oct12-2012.xlsx",sheet=3) %>%
  left_join(.,files1,by=c("SampleID"="file"))
#write_xlsx(covar,"~/Documents/Aldrin/CPDR_RNA-seq/Cell/covar.xlsx")


###Celltype
ct=read_csv("CIBERSORT.Output_Job12.csv")[,-1] %>% bind_cols(covar[,1],.) %>%column_to_rownames('SampleID')
ct=ct[,1:22]



ct=read_csv("CIBERSORT.Output_Abs_Job12.csv")[,-1] %>% bind_cols(covar[,1],.) %>%column_to_rownames('SampleID')
ct=ct[,1:22]

rn_covar<-covar %>%column_to_rownames("SampleID") %>% select(`Sample Type`,Race,`ERG Status`)





##Try ward.D, ward.D2, complete, ceneroid
##EUCLIDEAN

pheat1=pheatmap(ct,annotation_row = rn_covar,
                show_rownames = F,show_colnames=T,
                clustering_method="ward.D",
                clustering_distance_rows  = "euclidean")

pheat1


##Read in 

est_matrix=read_csv("../Cell/estimation_matrix.csv")
names(est_matrix)[-1]=covar$`FP #`
names(est_matrix)[1]="CellType"


covar_race=covar %>% slice(1:14)
x=est_matrix[1,c(-1,-15:-18)]
s=apply(est_matrix[,c(-1,-15:-18)],1,function(x){
  summary(lm(unlist(x)~covar_race$`ERG Status`))$coefficients[2,4]
})
which(s<.1)
unlist(est_matrix[,1])[which(s<.1)]
est_matrix %>% select(1:15)
