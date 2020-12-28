norm=["/home/dnousome/Documents/Armstrong/MD/QB/SH5519_SA35822_S89_L003.md.bam",'/home/dnousome/Documents/Armstrong/MD/QB/SH5526_SA35823_S112_L003.md.bam']
tum=["/home/dnousome/Documents/Armstrong/MD/TAGC/43-07647.md.bam", "/home/dnousome/Documents/Armstrong/MD/TAGC/43-07637.md.bam"]
normname=["SH5519_SA35822_S89_L003", "SH5526_SA35823_S112_L003"]
out= ["SH5519_SA35822_S89_L003.som.vcf.gz", "SH5526_SA35823_S112_L003.som.vcf.gz"]


configfile: "config.yaml"
localrules: all
rule target:
    input: 
       out


rule mutect_matched:
    input:
        tumor=tum,
        normal=norm
    output:
        "{out}"
    conda:
        'gatk4.yml'
    params:
        hgfa = config['hgfa'],
        illuminaInt= config['illuminaInt'],
        gnomadvcf= config['gnomadvcf'],
        samplename= normname
    shell:
        """
        gatk Mutect2 -R {params.hgfa} -I {input.tumor} -I {input.normal} -normal {params.samplename} --germline-resource {params.gnomadvcf} -O {output} --tmp-dir ~/Documents/Armstrong/temp
        """
