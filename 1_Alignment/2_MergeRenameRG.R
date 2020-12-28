####################
#
#Merge all the QB/TAGC multiples together
#
####################


##4 duplicates


qb1<-Sys.glob("~/Documents/Armstrong/Scratch/QB/*")
qbi<-gsub("FT-","",sapply(strsplit(basename(Sys.glob("~/Documents/Armstrong/Scratch/QB/*")),"_"),'[',2))
unique(qb)

qb1<-qb1[duplicated(qbi)|duplicated(qbi,fromLast=T)]
qbi<-qbi[duplicated(qbi)|duplicated(qbi,fromLast=T)]

qbs<-split(qb1,qbi)
code<-sapply(names(qbs),function(x){
  sprintf("samtools merge /home/dnousome/Documents/Armstrong/Scratch/QB/%s.prg.aligned.bam %s -@ 18",
  x,paste(qbs[[x]],collapse=" "))
  
})

code1<-c("#!/bin/bash",code)
write.table(code1,"Documents/Aldrin/APOLLO3/WES/Alignment/merge_qb.sh",row.names=F,quote=F,col.names=F)



#REPLACE RGs
code<-sapply(names(qbs),function(x){
  sprintf("gatk AddOrReplaceReadGroups -I /home/dnousome/Documents/Armstrong/Scratch/QB/%s.prg.aligned.bam -RGLB %s -RGPL %s -RGPU %s -RGSM %s -O /home/dnousome/Documents/Armstrong/Scratch/QB/%s.aligned.bam --TMP_DIR=/home/dnousome/Documents/Armstrong/temp",
          x,x,x,x,x,x)
  
})

code1<-c("#!/bin/bash",code)
write.table(code1,"Documents/Aldrin/APOLLO3/WES/Alignment/2_merge_qb_rg.sh",row.names=F,quote=F,col.names=F)





#########DO THE SAME FOR TAGC
tag<-Sys.glob("~/Documents/Armstrong/Scratch/TAGC/*")
tagi<-gsub("FT-","",sapply(strsplit(basename(Sys.glob("~/Documents/Armstrong/Scratch/TAGC//*")),"_"),'[',1))
unique(tagi)

tag<-tag[duplicated(tagi)|duplicated(tagi,fromLast=T)]
tagi<-tagi[duplicated(tagi)|duplicated(tagi,fromLast=T)]


tags<-split(tag,tagi)
code<-sapply(names(tags),function(x){
  sprintf("samtools merge /home/dnousome/Documents/Armstrong/Scratch/TAGC/%s.prg.aligned.bam %s -@ 18",
          x,paste(tags[[x]],collapse=" "))
  
})

code1<-c("#!/bin/bash",code)
write.table(code1,"Documents/Aldrin/APOLLO3/WES/Alignment/merge_tagc.sh",row.names=F,quote=F,col.names=F)






#REPLACE RGs
code<-sapply(names(tags),function(x){
  sprintf("gatk AddOrReplaceReadGroups -I /home/dnousome/Documents/Armstrong/Scratch/TAGC/%s.prg.aligned.bam -RGLB %s -RGPL %s -RGPU %s -RGSM %s -O /home/dnousome/Documents/Armstrong/Scratch/TAGC/%s.aligned.bam --TMP_DIR=/home/dnousome/Documents/Armstrong/temp",
          x,x,x,x,x,x)
  
})



code1<-c("#!/bin/bash",code)
write.table(code1,"Documents/Aldrin/APOLLO3/WES/Alignment/2_merge_tagc_rg.sh",row.names=F,quote=F,col.names=F)





###########

#16 unique only
##
java -jar picard.jar MergeSamFiles \
I=input_1.bam \
I=input_2.bam \
O=output_merged_files.bam

###TRY THIS 
samtools merge test_samtools.bam SH5518_SA35801_S41_L001.aligned.bam SH5518_SA35801_S41_L002.aligned.bam  -@ 18

conda activate gatk4
gatk AddOrReplaceReadGroups -I test_gatk.bam -RGLB 'SH5518_SA35801_S41_L001' -RGPL 'SH5518_SA35801_S41_L001' -RGPU 'SH5518_SA35801_S41_L001' -RGSM 'SH5518_SA35801_S41_L001' -O test_gatk1.bam --TMP_DIR=~/Documents/Armstrong/temp 
gatk AddOrReplaceReadGroups -I test_samtools.bam -RGLB 'SH5518_SA35801_S41_L001' -RGPL 'SH5518_SA35801_S41_L001' -RGPU 'SH5518_SA35801_S41_L001' -RGSM 'SH5518_SA35801_S41_L001' -O test_samtools1.bam --TMP_DIR=~/Documents/Armstrong/temp 
