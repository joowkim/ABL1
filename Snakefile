import pandas as pd
import numpy as np
import os


## load configure info
configfile: "configs/configs.yaml"


## load sample info
# samplesheet="./samples.tsv"
# units = pd.read_table(samplesheet, dtype={"sample" : str, "sample_group" : str })


def get_sample_name():
    import glob
    sample_list = glob.glob("00.rawdata/*.fq")
    sample_name = [os.path.basename(i).split(".fq")[0] for i in
                   sample_list]
    return sample_name


rule all:
    input:
        # expand("analysis/samtools/{sample}.flagstats.tsv",sample=get_sample_name()),
        "analysis/multiqc/multiqc_report.html",
        expand("analysis/bcftools/{sample}.vcf",sample=get_sample_name()),
        expand("analysis/freebayes/{sample}.vcf",sample=get_sample_name()),
        expand("analysis/get_coverage/{sample}.txt",sample=get_sample_name()),
        expand("analysis/get_coverage/{sample}.tsv",sample=get_sample_name()),
        expand("analysis/plot_coverage/{sample}.coverage.tsv",sample=get_sample_name()),
        expand("analysis/plot_coverage/{sample}.png",sample=get_sample_name()),
        expand("analysis/haplotype_caller/{sample}.vcf",sample=get_sample_name()),
        expand("analysis/mutect2/{sample}.vcf",sample=get_sample_name()),

# expand("analysis/09.samtools/{sample}.bam", sample=get_sample_name()),


# rule demux:
#     input:
#     output:
#     shell:


rule fastqc:
    """
    Run fastqc on raw_data/files.
    """
    input:
        "00.rawdata/{sample}.fq"
    output:
        html="analysis/fastqc/{sample}_fastqc.html",
        zip="analysis/fastqc/{sample}_fastqc.zip"
    params:
        outdir="analysis/fastqc/"
    log:
        stdout="logs/fastqc/{sample}.o",
        stderr="logs/fastqc/{sample}.e"
    threads: 1
    resources:
        mem_gb=8
    container:
        config["fastqc"]
        # "docker://biocontainers/fastqc:v0.11.9_cv8"
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """


rule fastq_screen:
    input:
        "00.rawdata/{sample}.fq",
    output:
        html="analysis/fastq_screen/{sample}_screen.html",
        txt="analysis/fastq_screen/{sample}_screen.txt",
    log:
        "logs/fastq_screen/{sample}.log",
    params:
        fastq_screen=config["fastq_screen"]["path"],
        conf=config["fastq_screen"]["conf"],
        bowtie2=config["bowtie2"]["path"],
    threads: 2
    container:
        config["fastq_screen"]["sif"]
        # "docker://quay.io/biocontainers/fastq-screen:0.14.0--pl5321hdfd78af_2"
    shell:
        """
        {params.fastq_screen} --aligner bowtie2 --bowtie2 {params.bowtie2} --outdir analysis/fastq_screen/ --threads {threads} --conf {params.conf}  {input} 2> {log}
        """


rule cutadapt:
    """
    Run cutadapt to remove short reads/adapter - less than 100bp 
    """
    input:
        "00.rawdata/{sample}.fq",
    output:
        trimmed_fq="analysis/cutadapt/{sample}.fq",
    log:
        stdout="logs/cutadapt/{sample}.o",
        stderr="logs/cutadapt/{sample}.e",
    params:
        length="200",
        adapter="GTAC",
    threads: 4
    resources:
        mem_gb=8
    container:
        config["cutadapt"],
        # "docker://quay.io/biocontainers/cutadapt:4.1--py38hbff2b2d_1"
    shell:
        """
        cutadapt -o {output.trimmed_fq}  --minimum-length {params.length}  -b {params.adapter}  --revcomp {input} > {log.stdout} 2> {log.stderr}
        """


rule pollux:
    """
    Run pollux to correct homopolymer errors.
    """
    input:
        rules.cutadapt.output.trimmed_fq,
    output:
        fqgz="analysis/pollux/{sample}.fq.gz",
        tempfq=temp("analysis/pollux/{sample}.fq.corrected"),
        lowfq=temp("analysis/pollux/{sample}.fq.low"),
    params:
        outdir="analysis/pollux",
    shell:
        """
        pollux -i {input} -o {params.outdir} -h "true"
        
        gzip -c {output.tempfq} > {output.fqgz}
        
        """


rule bwa:
    """
    Run bwa-mem
    """
    input:
        fq=rules.pollux.output.fqgz,
        idx=config['ref_modi']['index'],
    output:
        outsam=temp("analysis/bwamem/{sample}.sam")
    params:
        prefix="{sample}",
        # idx=config['ref_modi']['index'],
    log:
        stdout="logs/bwamem/{sample}.o",
        stderr="logs/bwamem/{sample}.e",
    threads: 4
    resources:
        mem_gb=20
    container:
        config["bwa"],
        # "docker://biocontainers/bwa:v0.7.17_cv1"
    shell:
        """
        bwa mem -M -t {threads} -R '@RG\\tID:{params.prefix}\\tLB:{params.prefix}\\tPL:Ion\\tPM:Torren\\tSM:{params.prefix}' {input.idx} {input.fq} > {output.outsam}  
        """

rule sam2bam:
    input:
        insam=rules.bwa.output.outsam,
    output:
        outbam="analysis/bwamem/{sample}.bam",
        outbai="analysis/bwamem/{sample}.bam.bai",
        idxstat="analysis/bwamem/{sample}.bam.idxstat",
    log:
        stdout="logs/sam2bam/{sample}.o",
        stderr="logs/sam2bam/{sample}.e",
    threads: 4
    resources:
        mem_gb=20
    container:
        config["samtools"],
    shell:
        """
        samtools view -S -b {input.insam} | samtools sort -O BAM -o {output.outbam} 1> {log.stdout} 2> {log.stderr}
        samtools index -@ {threads} {output.outbam}
        samtools idxstats {output.outbam} > {output.idxstat}
        """



rule samtools_stat:
    input:
        rules.sam2bam.output.outbam,
    output:
        out_stats="analysis/samtools_stats/{sample}.stats",
        out_flagstats="analysis/samtools_stats/{sample}.flagstat.txt",
    threads: 2
    container:
        config["samtools"],
        # "docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    shell:
        """
        samtools stats {input} > {output.out_stats}
        samtools flagstat -O tsv {input} > {output.out_flagstats}
        """


rule samtools:
    """
    Filter out supplementary alignment reads
    """
    input:
        bam=rules.sam2bam.output.outbam,
        outbai=rules.sam2bam.output.outbai,
        idxstat=rules.sam2bam.output.idxstat,
    output:
        output_coverage_tsv="analysis/plot_coverage/{sample}.coverage.tsv",
        tempbam=temp("analysis/samtools{sample}.temp.bam"),
        outbam=("analysis/samtools/{sample}.bam"),
        idxstat="analysis/samtools/{sample}.bam.idxstat",
        bam_idx="analysis/samtools/{sample}.bam.bai",
    log:
        stdout="logs/samtools/{sample}.o",
        stderr="logs/samtools/{sample}.e",
    # params:
    #     temp_bam=temp("analysis/03.samtools/SA.only.bam"),
    #     temp_id_list=temp("analysis/03.samtools/SA.id.list"),
    threads: 4
    resources:
        mem_gb=4
    container:
        config["samtools"],
        # "docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    shell:
        """
        samtools view --expr '![SA]' -O BAM -o {output.tempbam}  {input.bam}
        
        # If you want to keep uniquely-mapped reads + MAPQ > 20
        # -F 2308
        # read unmapped, not primary and supplementary
        samtools sort -O BAM {output.tempbam} | samtools view -h -F 2308 -q 20 -o {output.outbam} -O BAM
    
        samtools index {output.outbam}
        samtools idxstats {output.outbam} > {output.idxstat}
        samtools depth {output.outbam} > {output.output_coverage_tsv}
        """

rule get_align_metrics:
    """
    Run picard to get alignment metrics
    https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering
    """
    input:
        inbam=rules.samtools.output.outbam
    output:
        metrics="analysis/align_metrics/{sample}.align.metrics.txt",
    params:
        picard=config["picard"],
        ref_fa=config["ref_modi"]["sequence"],
    # container:
    #     "docker://quay.io/biocontainers/picard:2.26.10--hdfd78af_0"
    shell:
        """
        java -jar {params.picard} CollectAlignmentSummaryMetrics         I={input.inbam}         R={params.ref_fa}          O={output.metrics}         
        #METRIC_ACCUMULATION_LEVEL=READ_GROUP \
        #METRIC_ACCUMULATION_LEVEL=SAMPLE \
        """


rule mark_dups:
    """
    Run picard to mark duplicated reads
    https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering
    """
    input:
        inbam=rules.samtools.output.outbam,
        metrics=rules.get_align_metrics.output.metrics,
    output:
        outbam="analysis/mark_dups/{sample}.bam",
        outbamidx="analysis/mark_dups/{sample}.bam.bai",
    log:
        out="logs/mark_dup/{sample}.o",
        err="logs/mark_dup/{sample}.e"
    params:
        picard=config["picard"],
        ref_fa=config["ref_modi"]["sequence"],
        tmp_dir="/tmp",
    shell:
        #TAGGING_POLICY=All
        """
        java -jar {params.picard} MarkDuplicates         TMP_DIR={params.tmp_dir}         I={input.inbam}         O={output.outbam}         METRICS_FILE={input.metrics}         REMOVE_DUPLICATES=false 1>{log.out} 2>{log.err}
        samtools index {output.outbam}
        """


rule gatk_haplotype_caller:
    """
    Run gatk haplotype caller
    """
    input:
        inbam=rules.samtools.output.outbam,
    output:
        outvcf="analysis/haplotype_caller/{sample}.vcf",
    params:
        #option="-Xmx4g -XX:ParallelGCThreads=1",
        ref_fa=config["ref_modi"]["index"],
    # container:
    #     "docker://quay.io/biocontainers/gatk4:4.2.4.0--hdfd78af_0"
    shell:
        """
        gatk HaplotypeCaller -R {params.ref_fa} -I {input.inbam} -O {output.outvcf}
        """


rule gatk_mutect2:
    """
    Run gatk mutect2
    """
    input:
        inbam=rules.samtools.output.outbam,
    output:
        outvcf="analysis/mutect2/{sample}.vcf",
    params:
        ref_fa=config["ref_modi"]["index"],
    # container:
    #     "docker://quay.io/biocontainers/gatk4:4.2.4.0--hdfd78af_0"
    shell:
        """
          gatk Mutect2 -R {params.ref_fa} -I {input.inbam}    -O {output.outvcf}
        """


rule multiqc:
    """
    Run multiqc
    """
    input:
        # expand("analysis/star/{sample}.Log.final.out", sample=get_sample_name()),
        expand("analysis/fastqc/{sample}_fastqc.html",sample=get_sample_name()),
        expand("analysis/fastqc/{sample}_fastqc.zip",sample=get_sample_name()),
        # expand("analysis/bwamem/{sample}.bam.idxstat", sample=get_sample_name()),
        # expand("analysis/samtools/{sample}.bam.idxstat",sample=get_sample_name()),
        # expand("analysis/fastp/{sample}.fastp.json",sample=get_sample_name()),
        # expand("analysis/samtools/{sample}.flagstats.tsv",sample=get_sample_name()),
        # expand("analysis/align_metrics/{sample}.align.metrics.txt",sample=get_sample_name()),
        # expand("logs/mark_dup/{sample}.o",sample=get_sample_name()),
        expand("analysis/samtools_stats/{sample}.stats",sample=get_sample_name()),
        # expand("analysis/fastq_screen/{sample}_screen.txt",sample=get_sample_name()),
        # expand("analysis/samtools_stats/{sample}.flagstat.txt", sample=get_sample_name()),
    output:
        "analysis/multiqc/multiqc_report.html",
    log:
        stdout="logs/multiqc/multiqc.o",
        stderr="logs/multiqc/multiqc.e",
    threads: 4
    resources:
        mem_gb=100
    container:
        config["multiqc"],
        # "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0"
    shell:
        """
        multiqc -f {input} \
        -o analysis/multiqc \
        -n multiqc_report.html \
        --cl-config 'max_table_rows: 999999' \
        """


rule bcftools:
    """
    Run bcftools to call variants
    """
    input:
        inbam=rules.samtools.output.outbam
    output:
        outvcf="analysis/bcftools/{sample}.vcf"
    params:
        ref_fa=config['ref_modi']['index'],
    container:
        config['bcftools'],
        # "docker://quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0"
    shell:
        """
        bcftools mpileup -f {params.ref_fa} {input.inbam} | bcftools call -mv -Ov -o {output.outvcf}
        """


rule freebayes:
    """
    Run freebayes to call variants
    """
    input:
        inbam=rules.samtools.output.outbam
    output:
        outvcf="analysis/freebayes/{sample}.vcf"
    params:
        ref_fa=config['ref_modi']['index'],
    container:
        config["freebayes"]
        # "docker://quay.io/biocontainers/freebayes:1.3.6--hb089aa1_0"
    shell:
        """
        freebayes -f {params.ref_fa} -C 10 {input.inbam} > {output.outvcf}
        """


rule get_coverage:
    """
    Run David's script
    """
    input:
        inbam=rules.samtools.output.outbam,
    output:
        summaryfile="analysis/get_coverage/{sample}.txt",
        bedfile="analysis/get_coverage/{sample}.tsv",
    shell:
        """
        python3 src/get_coverage.py -b {input.inbam} -o {output.bedfile}  -min_read_over_lap 0.5  -min_amplicon_cov 300  -min_sample_cov 100  -min_read_uniform 0.2 -s {output.summaryfile}
        """


rule plot_coverage:
    """
    Run coverage.R
    """
    input:
        inbam=rules.samtools.output.outbam,
        coverage_tsv=rules.samtools.output.output_coverage_tsv,
    output:
        # output_coverage_tsv="analysis/plot_coverage/{sample}.coverage.tsv",
        outpng="analysis/plot_coverage/{sample}.png",
    params:
        prefix="analysis/plot_coverage/{sample}"
    container:
        config["tidyverse"]
    shell:
        """
        Rscript src/coverage.R {input.coverage_tsv} {params.prefix}
        """


#deprecated. after filtering, many reads filtered out.
rule fastp:
    """
    Run fastp - this is a quality trimming step.
    """
    input:
        "00.rawdata/{sample}.fq"
    output:
        trimmed_fq_gz="analysis/fastp/{sample}.fq.gz",
        outjson="analysis/fastp/{sample}.fastp.json",
        outhtml="analysis/fastp/{sample}.fastp.html",
    log:
        stdout="logs/fastp/{sample}.o",
        stderr="logs/fastp/{sample}.e",
    threads: 4
    resources:
        mem_gb=8
    container:
        config["fastp"],
        # "docker://quay.io/biocontainers/fastp:0.23.2--h79da9fb_0"
    shell:
        """
        fastp -i {input} -o {output.trimmed_fq_gz} --json {output.outjson} --html {output.outhtml} --thread 3 --adapter_sequence GTAC
        """


rule STAR:
    input:
        rules.cutadapt.output.trimmed_fq
    output:
        # see STAR manual for additional output files
        bam="analysis/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="analysis/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        log_final="analysis/star/{sample}.Log.final.out",
        log="analysis/star/{sample}.Log.out",
        rpg="analysis/star/{sample}.ReadsPerGene.out.tab",
        sj="analysis/star/{sample}.SJ.out.tab",
        g_dir=directory("analysis/star/{sample}._STARgenome"),
        pass1_dir=directory("analysis/star/{sample}._STARpass1"),
    params:
        # path to STAR reference genome index
        index=config["ref_modi"]["star_index"],
        outprefix="analysis/star/{sample}."
    log:
        "logs/star/{sample}.log"
    threads: 8
    resources:
        nodes=1,
        mem_gb=50,
    container:
        config["star"],
        # "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input} \
        --twopassMode Basic \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.outprefix} \
        --quantMode GeneCounts \
        --outStd Log 2> {log}
        samtools index {output.bam}
        """


rule star_samtools:
    input:
        rules.STAR.output.bam,
    output:
        outbam="analysis/star_samtools/{sample}.bam",
        bamidx="analysis/star_samtools/{sample}.bam.bai"
    threads: 4
    log:
        stdout="logs/star_samtools/{sample}.o",
        stderr="logs/star_samtools/{sample}.e",
    container:
        "docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    shell:
        """
        # If you want to keep uniquely-mapped reads + MAPQ > 20 #
        samtools view -h -F 256 -q 20 {input} | samtools sort -O BAM -o {output.outbam} > {log.stdout} 2> {log.stderr}
        
        samtools index -@ {threads} {output.outbam} 
        """
