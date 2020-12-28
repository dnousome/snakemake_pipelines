########################
#
#
#CNV Analysis
#
#
########################

##Step 0
gatk PreprocessIntervals \
-L nextera_hg38_1.bed \
-R ~/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
-O cnv_hg38.interval_list



############Create the PANEL OF NORMALS with all 28 samples
setwd("~/Documents/Armstrong/4_CNV")
norms=Sys.glob("*norm.counts.hdf5")
####Remove the duplicate norm samples
norms=norms[-22]
normscat=paste(paste0("-I ",norms),collapse=" ")
normout=sprintf("gatk --java-options -Xmx12g CreateReadCountPanelOfNormals %s --minimum-interval-median-percentile 5.0 -O cnvponC.pon.hdf5",normscat)

normout=c("#!/bin/bash",normout)
write.table(normout,"normout.sh",row.names=F,col.names=F,quote=F)
