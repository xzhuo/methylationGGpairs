import os

SUFFIX = '.bam'
suffix_length = len(SUFFIX)
SAMPLES = set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir("."))))
# SAMPLES = ["bam"]
print(SAMPLES)

REF = "/storage1/fs1/hprc/Active/xzhuo/ref/hg38.fa"
WGBS = "/storage1/fs1/hprc/Active/xzhuo/pacbio_ggpairs/HG002.wgbs.methylC.CpG.bed"

rule all:
    input:
        expand("ggpairs/HG002.wgbs.{sample}.ggpairs.pdf", sample=SAMPLES),
        expand("ggpairs/HG002.wgbs.{sample}.ggpairs.txt", sample=SAMPLES)

rule pbmm2:
    input:
        bam = "ubam/{sample}.bam",
        ref = REF
    output:
        bam = "hg38/{sample}.hg38.bam",
        bai = "hg38/{sample}.hg38.bam.bai"
    threads:
        8
    # container:
    #     "docker://xiaoyuz/biotools:latest"
    shell:
        "mkdir -p hg38;pbmm2 align {input.ref} {input.bam} {output.bam} --preset CCS --sort -j {threads} -m 2G"

rule pbCpGtools:
    input:
        bam = "hg38/{sample}.hg38.bam"
    output:
        model = "pbcpgtools/{sample}.model.combined.bed",
        counts = "pbcpgtools/{sample}.count.combined.bed"
    threads:
        8
    params:
        model = "/opt/conda/pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
        prefix_model = "pbcpgtools/{sample}.model",
        prefix_count = "pbcpgtools/{sample}.count"
    # container:
    #     "docker://xiaoyuz/modbamutil:latest"
    shell:
        "mkdir -p pbcpgtools;aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_model} --model {params.model} --threads {threads}; \
            aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_count} --pileup-mode count --threads {threads}"

rule ggpairs_wgbs:
    input:
        model = "pbcpgtools/{sample}.model.combined.bed",
        counts = "pbcpgtools/{sample}.count.combined.bed",
        wgbs = WGBS
    output:
        bed = "ggpairs/HG002.wgbs.{sample}.methylC.CpG.bed",
        txt = "ggpairs/HG002.wgbs.{sample}.ggpairs.txt"
    threads:
        4
    # container:
    #     "docker://xiaoyuz/biotools:latest"
    shell:
        """mkdir -p ggpairs;bedtools intersect -a {input.wgbs} -b {input.counts} -loj|cut -f1-5,9,11 | \
            bedtools intersect -a stdin -b {input.model} -loj |cut -f1-7,11,13 > {output.bed}; \
            perl -lane 'print join("\\t",$F[0],$F[1],$F[2],$F[3],$F[5],$F[7]) if $F[4]>=5 && $F[6] >=5 && $F[8]>=5' {output.bed} > {output.txt}"""

rule ggpairs:
    input:
        "ggpairs/HG002.wgbs.{sample}.ggpairs.txt"
    output:
        "ggpairs/HG002.wgbs.{sample}.ggpairs.pdf"
    # conda:
    #     "envs/ggplot.yaml"
    threads:
        4
    script:
        "scripts/methylation.ggpairs.R"