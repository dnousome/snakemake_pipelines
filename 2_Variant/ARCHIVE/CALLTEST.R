##create final list
se

#
##No INT 9;50
##55 mins40s
time gatk Mutect2 -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-I SH5529_SA35828_S122_L003.brecal.bam -max-mnp-distance 0 -O normal_test_noint.vcf.gz --tmp-dir ~/Documents/Armstrong/temp
time gatk Mutect2 -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna -I SH5529_SA35824_S121_L003.brecal.bam -max-mnp-distance 0 -O normal_test1_noint.vcf.gz --tmp-dir ~/Documents/Armstrong/temp

gatk CombineGVCFs -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-L /home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed \
--variant normal_test_noint_header.vcf.gz \
--variant normal_test1_noint_header.vcf.gz \
-O cohort.vcf.gz

##INFO=<ID=END,Number=1,Type=Integer,Description="DUMMY">

bcftools view -h normal_test_noint.vcf.gz|tail -n1 >samp.header
cat newheader1 header2 samp.header >newheader2

bcftools reheader normal_test_noint.vcf.gz -h newheader2 -o normal_test_noint_header.vcf.gz


bcftools view -h normal_test1_noint.vcf.gz|tail -n1 >samp.header
cat newheader1 header2 samp.header >newheader2


bcftools reheader normal_test1_noint.vcf.gz -h newheader2 -o normal_test1_noint_header.vcf.gz




gatk CreateSomaticPanelOfNormals -V cohort.vcf.gz \
-O pon_noint.vcf.gz





###INT FILE
##INT? #1050

time gatk Mutect2 -L /home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna -I SH5529_SA35828_S122_L003.brecal.bam -max-mnp-distance 0 -O normal_test_int.vcf.gz --tmp-dir ~/Documents/Armstrong/temp
time gatk Mutect2 -L /home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna -I SH5529_SA35824_S121_L003.brecal.bam -max-mnp-distance 0 -O normal_test1_int.vcf.gz --tmp-dir ~/Documents/Armstrong/temp


bcftools view -h normal_test_int.vcf.gz|tail -n1 >samp.header
cat newheader1 header2 samp.header >newheader2
bcftools reheader normal_test_int.vcf.gz -h newheader2 -o normal_test_int_header.vcf.gz


bcftools view -h normal_test1_int.vcf.gz|tail -n1 >samp.header
cat newheader1 header2 samp.header >newheader2
bcftools reheader normal_test1_int.vcf.gz -h newheader2 -o normal_test1_int_header.vcf.gz

for i in *_header.vcf.gz; do
bcftools index $i -t; done
gatk CombineGVCFs -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-L /home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed \
--variant normal_test_int_header.vcf.gz \
--variant normal_test1_int_header.vcf.gz \
-O cohort1.vcf.gz


gatk CreateSomaticPanelOfNormals -V cohort1.vcf.gz \
-O pon_int.vcf.gz






gatk Mutect2 \
-R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-I SA35801.brecal.bam \
--germline-resource ~/dn/Annotations/vcfs/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
--panel-of-normals ~/Documents/Aldrin/APOLLO3/WES/Alignment/pon.vcf.gz \
-O test_call.vcf.gz

gatk Mutect2 \
-R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-I SA35801.brecal.bam \
-I SA35815.brecal.bam \
-normal SA35815 \
--germline-resource ~/dn/Annotations/vcfs/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
-O test_matched.vcf.gz

###Compare the matching?
gatk Mutect2 \
-R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-I TAGC/43-07647.brecal.bam \
--germline-resource ~/dn/Annotations/vcfs/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
--panel-of-normals ~/Documents/Aldrin/APOLLO3/WES/Alignment/pon.vcf.gz \
-O test_call_TAGConly.vcf.gz

gatk Mutect2 \
-R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-I TAGC/43-07647.brecal.bam \
-I QB/SH5519_SA35822_S89_L003.brecal.bam \
-normal SH5519_SA35822_S89_L003 \
--germline-resource ~/dn/Annotations/vcfs/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
-O test_matched_QBTAGC.vcf.gz




snakemake -s mutect2.snake --config tumor="/home/dnousome/Documents/Armstrong/RECAL/TAGC/43-07647.brecal.bam"  output="NEW.vcf" normal="/home/dnousome/Documents/Armstrong/RECAL/QB/SH5519_SA35822_S89_L003.brecal.bam" samplename="SA35822" --use-conda -np 



gatk Mutect2 -L "/home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed" \
--germline-resource /home/dnousome/dn/Annotations/vcfs/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
-R "/home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna" \
-I ../RECAL/TAGC/43-07647.brecal.bam -max-mnp-distance 0 -O 43-07647.test_GERM.vcf.gz \
--tmp-dir ~/Documents/Armstrong/temp

gatk Mutect2 -L "/home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed" \
--germline-resource /home/dnousome/dn/Annotations/vcfs/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
-R "/home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna" \
-I ../RECAL/QB/SH5519_SA35822_S89_L003.brecal.bam -max-mnp-distance 0 -O SA35822.test_GERM.vcf.gz \
--tmp-dir ~/Documents/Armstrong/temp


FT-SA35822	
43-07647	


##Cvoerage
bedtools coverage -a ~/dn/Annotations/illumina/nextera_hg38_1.bed -b SH5528_SA35793_S117_L003.brecal.bam -mean >testCOVER

awk '{total+=$4} END {print total/NR}' testCOVER 

gatk CreateSequenceDictionaryTool
gatk BedToIntervalList -I "/home/dnousome/dn/Annotations/illumina/nextera_hg38_1.bed" \
-O "/home/dnousome/dn/Annotations/illumina/nextera_hg38_1.interval" \
-SD /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.dict
gatk CollectHsMetrics -I  SH5528_SA35793_S117_L003.brecal.bam -O test_GATKcover \
-TI ~/dn/Annotations/illumina/nextera_hg38_1.interval \
-BI ~/dn/Annotations/illumina/nextera_hg38_1.interval 
gatk CollectWgsMetrics -I  -O test_GATKcover -R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna

