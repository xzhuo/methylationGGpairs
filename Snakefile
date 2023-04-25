import os

SUFFIX = '.bam'
suffix_length = len(SUFFIX)
SAMPLES = set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir("."))))
# SAMPLES = ["bam"]
print(SAMPLES)

rule all:
    input:
        expand("methylation.{sample}.ggpairs.pdf", sample=SAMPLES),

rule pbmm2:
    input:
        bam = "{sample}.bam",
        ref = "/storage1/fs1/hprc/Active/xzhuo/ref/hg38.fa"
    output:
        bam = "{sample}.hg38.bam",
        bai = "{sample}.hg38.bam.bai"
    threads:
        8
    params:
        sort_threads = 2,
        alignment_threads = 6
    # container:
    #     "docker://xiaoyuz/biotools:latest"
    shell:
        "pbmm2 align {input.ref} {input.bam} {output.bam} --preset CCS --sort -j {params.alignment_threads} -J {params.sort_threads}"

rule pbCpGtools:
    input:
        bam = "{sample}.hg38.bam"
    output:
        model = "{sample}.model.combined.bed",
        counts = "{sample}.count.combined.bed"
    threads:
        8
    params:
        model = "/pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite",
        prefix_model = "{sample}.model",
        prefix_count = "{sample}.count"
    # container:
    #     "docker://xiaoyuz/modbamutil:latest"
    shell:
        "aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_model} --model {params.model} --threads {threads}; \
            aligned_bam_to_cpg_scores --bam {input} --output-prefix {params.prefix_count} --pileup-mode count --threads {threads}"

rule ggpairs_wgbs:
    input:
        model = "{sample}.model.combined.bed",
        counts = "{sample}.count.combined.bed",
        wgbs = "HG002.wgbs.methylC.CpG.bed"
    output:
        bed = "HG002.wgbs.{sample}.methylC.CpG.bed",
        txt = "HG002.wgbs.{sample}.ggpairs.txt"
    # container:
    #     "docker://xiaoyuz/biotools:latest"
    shell:
        """bedtools intersect -a {input.wgbs} -b {input.counts} -loj|cut -f1-5,9,10 | \
            bedtools intersect -a stdin -b {input.model} -loj |cut -f1-7,11,12 > {output.bed}"; \
            perl -lane 'print join("\t",$F[0],$F[1],$F[2],$F[3],$F[5],$F[7]) if $F[4]>=5 && $F[6] >=5 && $F[8]>=5' {output.bed} > {output.txt}"""

rule ggpairs:
    input:
        "HG002.wgbs.{sample}.{ref}.ggpairs.txt"
    output:
        "methylation.{sample}.{ref}.ggpairs.pdf"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/methylation.ggpairs.r"