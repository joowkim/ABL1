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
        # expand("analysis/04.samtools/{sample}.flagstats.tsv",sample=get_sample_name()),
        "analysis/06.multiqc/multiqc_report.html",
        expand("analysis/07.bcftools/{sample}.vcf",sample=get_sample_name()),
        expand("analysis/08.freebayes/{sample}.vcf",sample=get_sample_name()),
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
        html="analysis/01.fastqc/{sample}_fastqc.html",
        zip="analysis/01.fastqc/{sample}_fastqc.zip"
    params:
        outdir="analysis/01.fastqc/"
    log:
        stdout="logs/01.fastqc/{sample}.o",
        stderr="logs/01.fastqc/{sample}.e"
    threads: 1
    resources:
        mem_gb=8
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """


rule cutadapt:
    """
    Run cutadapt to remove short reads/adapter - less than 100bp 
    """
    input:
        "00.rawdata/{sample}.fq",
    output:
        trimmed_fq="analysis/02.cutadapt/{sample}.fq",
    log:
        stdout="logs/02.cutadapt/{sample}.o",
        stderr="logs/02.cutadapt/{sample}.e",
    params:
        length="100",
        adapter="GTAC",
    threads: 4
    resources:
        mem_gb=8
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
        fqgz="analysis/03.pollux/{sample}.fq.gz",
        tempfq=temp("analysis/03.pollux/{sample}.fq.corrected"),
        lowfq=temp("analysis/03.pollux/{sample}.fq.low"),
    params:
        outdir="analysis/03.pollux",
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
    output:
        outbam="analysis/04.bwamem/{sample}.bam",
        outbai="analysis/04.bwamem/{sample}.bam.bai",
        idxstat="analysis/04.bwamem/{sample}.bam.idxstat",
    params:
        prefix="{sample}",
        idx=config['ref_modi']['index'],
    log:
        stdout="logs/04.bwamem/{sample}.o",
        stderr="logs/04.bwamem/{sample}.e",
    threads: 4
    resources:
        mem_gb=20
    shell:
        """
        bwa mem -R '@RG\\tID:{params.prefix}\\tLB:{params.prefix}\\tPL:Ion\\tPM:Torren\\tSM:{params.prefix}' {params.idx} {input.fq}  | samtools sort -O BAM -o {output.outbam} > {log.stdout} 2> {log.stderr}
        
        samtools index -@ {threads} {output.outbam}
        samtools idxstats {output.outbam} > {output.idxstat}
        """


rule samtools:
    """
    Filter out supplementary alignment reads
    """
    input:
        bam=rules.bwa.output.outbam,
        outbai=rules.bwa.output.outbai,
        idxstat=rules.bwa.output.idxstat,
    output:
        tempbam=temp("analysis/05.samtools{sample}.temp.bam"),
        outbam=("analysis/05.samtools/{sample}.bam"),
        idxstat="analysis/05.samtools/{sample}.bam.idxstat",
        bam_idx="analysis/05.samtools/{sample}.bam.bai",
    log:
        stdout="logs/05.samtools/{sample}.o",
        stderr="logs/05.samtools/{sample}.e",
    # params:
    #     temp_bam=temp("analysis/03.samtools/SA.only.bam"),
    #     temp_id_list=temp("analysis/03.samtools/SA.id.list"),
    threads: 4
    resources:
        mem_gb=4
    shell:
        """
               
        samtools view --expr '![SA]' -O BAM -o {output.tempbam}  {input.bam}
        
        # If you want to keep uniquely-mapped reads + MAPQ > 20
        samtools sort -O BAM {output.tempbam} | samtools view -h -F 256 -q 20 -o {output.outbam} -O BAM
    

        samtools index {output.outbam}
        samtools idxstats {output.outbam} > {output.idxstat}
        """

rule multiqc:
    """
    Run multiqc
    """
    input:
        # expand("analysis/08.star/{sample}.Log.final.out", sample=get_sample_name()),
        expand("analysis/01.fastqc/{sample}_fastqc.html",sample=get_sample_name()),
        expand("analysis/01.fastqc/{sample}_fastqc.zip",sample=get_sample_name()),
        # expand("analysis/02.bwamem/{sample}.bam.idxstat", sample=get_sample_name()),
        expand("analysis/05.samtools/{sample}.bam.idxstat",sample=get_sample_name()),
    # expand("analysis/02.fastp/{sample}.fastp.json",sample=get_sample_name()),
    # expand("analysis/03.samtools/{sample}.flagstats.tsv",sample=get_sample_name()),
    output:
        "analysis/06.multiqc/multiqc_report.html",
    log:
        stdout="logs/06.multiqc/multiqc.o",
        stderr="logs/06.multiqc/multiqc.e",
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        multiqc -f {input} \
        -o analysis/06.multiqc \
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
        outvcf="analysis/07.bcftools/{sample}.vcf"
    params:
        ref_fa=config['ref_modi']['sequence'],
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
        outvcf="analysis/08.freebayes/{sample}.vcf"
    params:
        ref_fa=config['ref_modi']['sequence'],
    shell:
        """
        freebayes -f {params.ref_fa} -C 10 {input.inbam} > {output.outvcf}
        """


#deprecated. after filtering, no reads remained.
rule fastp:
    """
    Run fastp - this is a quality trimming step.
    """
    input:
        "00.rawdata/{sample}.fq"
    output:
        trimmed_fq_gz="analysis/02.fastp/{sample}.fq.gz",
        outjson="analysis/02.fastp/{sample}.fastp.json",
        outhtml="analysis/02.fastp/{sample}.fastp.html",
    log:
        stdout="logs/02.fastp/{sample}.o",
        stderr="logs/02.fastp/{sample}.e",
    threads: 4
    resources:
        mem_gb=8
    shell:
        """
        fastp -i {input} -o {output.trimmed_fq_gz} --json {output.outjson} --html {output.outhtml} --thread 3 --adapter_sequence GTAC
        """


rule STAR:
    input:
        rules.cutadapt.output.trimmed_fq
    output:
        # see STAR manual for additional output files
        bam="analysis/09.star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="analysis/09.star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        log_final="analysis/09.star/{sample}.Log.final.out",
        log="analysis/09.star/{sample}.Log.out",
        rpg="analysis/09.star/{sample}.ReadsPerGene.out.tab",
        sj="analysis/09.star/{sample}.SJ.out.tab",
        g_dir=directory("analysis/09.star/{sample}._STARgenome"),
        pass1_dir=directory("analysis/09.star/{sample}._STARpass1"),
    params:
        # path to STAR reference genome index
        index=config["ref_modi"]["star_index"],
        outprefix="analysis/09.star/{sample}."
    log:
        "logs/09.star/{sample}.log"
    threads: 8
    resources:
        nodes=1,
        mem_gb=50,
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
        outbam="analysis/10.samtools/{sample}.bam",
        bamidx="analysis/10.samtools/{sample}.bam.bai"
    threads: 4
    log:
        stdout="logs/10.samtools/{sample}.o",
        stderr="logs/10.samtools/{sample}.e",
    shell:
        """
        # If you want to keep uniquely-mapped reads + MAPQ > 20 #
        samtools view -h -F 256 -q 20 {input} | samtools sort -O BAM -o {output.outbam} > {log.stdout} 2> {log.stderr}
        
        samtools index -@ {threads} {output.outbam} 
        """

