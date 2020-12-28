

gatk Mutect2 \
-R /home/dnousome/dn/Annotations/ucsc_hg38/GCA_000001405.15_GRCh38_full_analysis_set.fna \
-I SA35801.brecal.bam \
-I SA35815.brecal.bam \
-normal SA35815 \
--germline-resource ~/dn/Annotations/vcfs/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
-O test_matched.vcf.gz
