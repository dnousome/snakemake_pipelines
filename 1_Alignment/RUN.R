#################################
#
#Additional R Code
#
#
#################################

##Check MD5 sums
#md5sum *.gz > CPDR.md5sum

setwd("~/Documents/APOLLO/APOLLO3/WES")

md5<-Sys.glob("*/*.md5sum")


cpdr<-seq(1,22,2)
og<-seq(2,22,2)

s<-mapply(function(x,y){
  print(x)
  x1<-read_delim(md5[x],col_names = F,delim=" ") %>%
    mutate(X2=gsub(" ", "",X2))
  y1<-read_delim(md5[y],col_names = F,delim=" ") %>%
    mutate(X2=basename(X2)) %>%
    filter(X2 %in% x1$X2)
  tibble(totalN_og=nrow(y1),totalN_cpdr=nrow(x1),
         matchmd5=sum(x1$X1==y1$X1))
  
},cpdr,og,SIMPLIFY=F)


##########Create script for MErging
f1<-Sys.glob("~/Documents/Aldrin/APOLLO3/WES/QB/*/*_R1_001.fastq.gz")

f2<-gsub("R1","R2",f1)

id<-sapply(strsplit(basename(f1),"_"),'[',2)
  

id=ifelse(duplicated(id)|duplicated(id,fromLast=T),paste(id,sapply(strsplit(basename(f1),"_"),'[',4),sep="_"),id)
  
codes<-sprintf("bwa mem -t 8 -R '@RG\\tID:%s\\tSM:%s' ~/dn/Annotations/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz %s %s | samtools sort -@8 -o %s_output.bam -",
        id,id,f1,f2,id)
setwd("/home/dnousome/Documents/Aldrin/APOLLO3/WES/Alignment/")

s<-c("#!/bin/bash",codes)
write.table(s,"run.sh",quote=F,row.names=F,col.names=F)
