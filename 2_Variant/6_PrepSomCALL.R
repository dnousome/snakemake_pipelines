####################################
#
#
#Prep Somatic calling dataset
#
#
####################################
library(tidyverse)
tagcsamps<-readxl::read_xlsx("~/Documents/Aldrin/APOLLO3/Data/10152019_APOLLO3_WGS WES RNASeq Prot_APOLLO IDs.xlsx")
qbsamps<-readxl::read_xlsx("~/Documents/Aldrin/APOLLO3/Data/QD597 - 42 WES DNA Sample-Sequencing QC Data.xlsx",skip = 1) %>%
  arrange(`Customer ID`) %>%
  filter(!is.na(`Customer ID`))


allmd1<-Sys.glob("~/Documents/Armstrong/RECAL/*/*.bam")
allmd2<-Sys.glob("~/Documents/Armstrong/MD/*/*.bam")

allvcf<-Sys.glob("~/Documents/Armstrong/PON/*/*rehead.vcf.gz")

allmdnames<-sapply(strsplit(basename(allmd2),"_"),function(x){
  x[grepl("SA",x) | grepl('[0-9]',substr(x,1,1))]

  })

allmdnames<-gsub("FT-","",sapply(strsplit(allmdnames,"[.]"),'[',1))

allsamps1<-tibble(allmd1,allmd2,allvcf,allmdnames)


qbsamps1<-qbsamps %>% 
  select(AID=`Accession ID`,GID=`Customer ID`) %>%
  mutate(AID1=gsub("FT-","",AID)) 

tagcsamps1<-tagcsamps %>% select(TID=`WES Seq ID`) %>%
  na.omit %>%
  mutate(AID=sapply(strsplit(TID,"-"),function(x)paste0(x[1],"-",x[2])),
         GID=sapply(strsplit(TID,"-"),function(x)paste0(x[3],"-",x[4])),
         AID1=AID) %>%
  select(AID,GID,AID1) 

allsamps2<-bind_rows(qbsamps1,tagcsamps1) %>%
  full_join(.,allsamps1,by=c("AID1"="allmdnames")) %>%
  mutate(GID1=sapply(strsplit(GID,"-"),'[',1))




normals<-allsamps2 %>% filter(grepl("N",GID))


####Create list to pool
##28 normals
write.table(normals$allvcf,"~/Documents/Aldrin/APOLLO3/WES/2_Variant/vcf.list",col.names = F,quote=F,row.names=F)


###Run the normals
gatk CombineGVCFs -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-L /home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed \
--variant vcf.list \
-O cohort_pon.vcf.gz

###PON1 is with the Intervals included during the calling phase
gatk CreateSomaticPanelOfNormals -V cohort_pon.vcf.gz -O pon.vcf.gz


##Rename the files to prep for somatic calling
x=split(allsamps2,allsamps2$GID1)[[1]]
te<-lapply(split(allsamps2,allsamps2$GID1),function(x){
  norm<-x$allmd1[grep("N",x$GID)]
  norm2<-x$allmd2[grep("N",x$GID)]
  tum<-x$allmd1[grep("T",x$GID)]
  tum2<-x$allmd2[grep("T",x$GID)]
  
  namenorm<-gsub(".md.bam","",basename(x$allmd2[1]))
  tibble(norm,norm2,tum,tum2,namenorm)
  })




s<-bind_rows(te) %>%
  mutate(namenorm=ifelse(duplicated(namenorm),paste0(namenorm,"-1"),namenorm))


mv.1=sprintf("mv %s /home/dnousome/Documents/Armstrong/FINALBAM/%s.norm.bam",s$norm,s$namenorm)
mv.2=sprintf("mv %s /home/dnousome/Documents/Armstrong/FINALBAM/%s.tum.bam",s$tum,s$namenorm)
mv.3=sprintf("mv %s /home/dnousome/Documents/Armstrong/FINALBAM/%s.norm.bai",gsub(".bam",".bai",s$norm),s$namenorm)
mv.4=sprintf("mv %s /home/dnousome/Documents/Armstrong/FINALBAM/%s.tum.bai",gsub(".bam",".bai",s$tum),s$namenorm)


mv.5=sprintf("mv %s /home/dnousome/Documents/Armstrong/BAM_MD/%s.norm.bam",s$norm2,s$namenorm)
mv.6=sprintf("mv %s /home/dnousome/Documents/Armstrong/BAM_MD/%s.tum.bam",s$tum2,s$namenorm)
mv.7=sprintf("mv %s /home/dnousome/Documents/Armstrong/BAM_MD/%s.norm.bai",paste0(s$norm2,".bai"),s$namenorm)
mv.8=sprintf("mv %s /home/dnousome/Documents/Armstrong/BAM_MD/%s.tum.bai",paste0(s$tum2,".bai"),s$namenorm)


#"mv /home/dnousome/Documents/Armstrong/MD/TAGC/43-07647.md.bam /home/dnousome/Documents/Armstrong/BAM_MD/SH5519_SA35822_S89_L003.tum.bam"           
#1187
# "mv /home/dnousome/Documents/Armstrong/MD/TAGC/43-07637.md.bam /home/dnousome/Documents/Armstrong/BAM_MD/SH5526_SA35823_S112_L003.tum.bam"                 
#1405
#"mv /home/dnousome/Documents/Armstrong/MD/TAGC/43-07642.md.bam /home/dnousome/Documents/Armstrong/BAM_MD/SH5529_SA35824_S121_L003.tum.bam"      
#1472

sapply(mv.1,system)
sapply(mv.2,system)
sapply(mv.3,system)
sapply(mv.4,system)


sapply(mv.5,system)
sapply(mv.6,system)
sapply(mv.7,system)
sapply(mv.8,system)


###For one duplicate sample create duplicate matched norm!
system("cp SH5527_FT-SA35807_S114_L003.norm.bam SH5527_FT-SA35807_S114_L003-1.norm.bam")
system("cp SH5527_FT-SA35807_S114_L003.norm.bai SH5527_FT-SA35807_S114_L003-1.norm.bai")

write_tsv(s,"~/Documents/Aldrin/APOLLO3/WES/2_Variant/rename_bams1.tsv")




















###Create 
scripts<-sapply(split(allsamps2,allsamps2$GID1),function(x){
  x1<-x %>% arrange(GID) 
  sprintf("snakemake -s mutect2.snake --config tumor=%s normal=%s samplename=%s m2matchedout=%s --cores=18 -j --use-conda",
          x1$allmd[2],
          x1$allmd[1],
          gsub(".md.bam","",basename(x1$allmd[1])),
          paste0(x1$AID[1],".somatic.vcf.gz"))
})
scripts<-sapply(split(allsamps2,allsamps2$GID1),function(x){
  x1<-x %>% arrange(GID) 
  sprintf("snakemake -s mutect2.snake --config tumor=%s normal=%s samplename=%s m2matchedout=%s --cores=18 -j --use-conda",
          x1$allmd[2],
          x1$allmd[1],
          gsub(".md.bam","",basename(x1$allmd[1])),
          paste0(x1$AID[1],".somatic.vcf.gz"))
})


scripts1<-c("#!/bin/bash",scripts)
write.table(scripts1,"~/Documents/Aldrin/APOLLO3/WES/Variant/mutect2.sh",row.names=F,quote=F,col.names=F)


#snakemake -s mutect2.snake --config tumor="/home/dnousome/Documents/Armstrong/RECAL/TAGC/43-07647.brecal.bam"  output="NEW.vcf" normal=/home/dnousome/Documents/Armstrong/RECAL/QB/SH5519_SA35822_S89_L003.brecal.bam samplename="SA35822" m2matchedout="TEST.vcf" mutect_tumoronly=True mutect_matched=False --use-conda -np -j6


paste(basename(s$norm),collapse="','")
paste(basename(s$tum),collapse="','")
paste(s$namenorm,collapse="','")
paste(ifelse(duplicated(s$namenorm),paste0(s$namenorm,".som_1.vcf.gz"),paste0(s$namenorm,".som.vcf.gz")),collapse="','")

dirname(s$norm)
dirname(s$tum)
write_tsv(allsamps2,"~/Documents/Aldrin/APOLLO3/WES/Variant/samples.file")


scripts1<-c("#!/bin/bash",scripts)
write.table(scripts1,"~/Documents/Aldrin/APOLLO3/WES/Variant/mutect2.sh",row.names=F,quote=F,col.names=F)

setwd("~/Documents/Aldrin/APOLLO3/WES/Variant/Mutect2/")
sapply(sprintf("cp %s %s.norm.bam",s$norm,s$namenorm)[1:2],system)
sapply(sprintf("cp %s %s.tum.bam",s$tum,s$namenorm)[1:2],system)

system