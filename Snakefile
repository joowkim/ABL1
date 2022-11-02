import glob
from typing import List

import pandas as pd
import os


## load configure info
configfile: "configs/configs.yaml"

## load sample info
samp_info = pd.read_csv(os.path.join("ubam",config["run_id"],config["sample_info_csv"]))


def get_sample_name() -> List[str]:
    # import glob
    # sample_list = glob.glob("00.rawdata/*.fq")
    # sample_name = [os.path.basename(i).split(".fq")[0] for i in
    #                sample_list]
    sample_name = samp_info["sample_id"].to_list()
    return sample_name


rule all:
    input:
        # expand("analysis/samtools/{sample}.flagstats.tsv",sample=get_sample_name()),
        expand("analysis/{run_id}/multiqc/multiqc_report.html",run_id=config["run_id"]),
        expand("analysis/{run_id}/bcftools/{sample}.vcf",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/freebayes/{sample}.vcf",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/get_coverage/{sample}.txt",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/get_coverage/{sample}.tsv",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/plot_coverage/{sample}.coverage.tsv",sample=get_sample_name(),run_id=config[
            "run_id"]),
        expand("analysis/{run_id}/plot_coverage/{sample}.png",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/deepvariant/{sample}.vcf",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/deepvariant/{sample}.gvcf",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/deepvariant/{sample}.visual_report.html",sample=get_sample_name(),run_id=config[
            "run_id"]),
        # expand("analysis/{run_id}/haplotype_caller/{sample}.vcf",sample=get_sample_name(),run_id=config["run_id"]),
        # expand("analysis/{run_id}/mutect2/{sample}.vcf",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/vcftools/{sample}.filt.recode.vcf",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/vcf_html/{sample}.vcf.html",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/fastq_scan/{sample}.summary.csv",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/fastq_scan/done",run_id=config["run_id"]),

    # expand("analysis/{run_id}/fastq_scan/fastq_stat.csv",run_id=config["run_id"]),


rule bam2fastq:
    input:
        ubam=os.path.join("ubam",config["run_id"],config["ubam"]),
    output:
        "analysis/{run_id}/bam2fastq/nomatch_rawlib.basecaller.fastq",
    container:
        config["samtools"],
    message: "Convert ubam into a fastq file."
    log:
        stderr="analysis/{run_id}/logs/bam2fastq/{run_id}.e"
    shell:
        """
        samtools bam2fq {input} > {output} 2>{log.stderr}
        """


rule process_samp_info:
    input:
        default_barcode=config["default_barcodes"],
        infq=rules.bam2fastq.output,
        samp_info_csv=os.path.join("ubam",config["run_id"],config["sample_info_csv"]),
    output:
        # barcode_info=os.path.join("analysis", config["run_id"], "barcode_info", "barcode_info.tsv"),
        barcode_info="analysis/{run_id}/barcode_info/barcode_info.tsv",
    message: "Generate a meta info file for demux."
    run:
        barcode_df: pd = pd.read_csv(input.default_barcode,sep="\t")
        # barcode_df looks like this.
        # id  seq style adapter
        # 420   CGAACATATTC    S5    GTAC
        # 421  TTGGACTTATTC    S5    GTAC
        # 422  CGAGGCAATGAC    S5    GTAC
        # 423  TCGAGATTAATC    S5    GTAC
        # 424   TTCGCCAACAC    S5    GTAC
        # 425  TCGGCACGAATC    S5    GTAC

        samp_info: pd = pd.read_csv(input.samp_info_csv)
        # samp_df looks like this.
        # row,sample_id,barcode
        # 1,6806006070_DL,420
        # 2,1_per_minor_qc,421
        # 3,690702015_DO,422
        # 4,1_per_major_qc,423

        samp_df: pd = barcode_df.merge(samp_info,left_on="id",right_on="barcode")
        final_df: pd = samp_df[["sample_id", "seq", "style", "adapter"]]
        final_df.columns = ["ID", "SEQ", "Style", "Adapter"]
        final_df.to_csv(output.barcode_info,index=False,sep="\t")


rule demux:
    input:
        barcode_info=rules.process_samp_info.output.barcode_info,
    output:
        touch=touch("analysis/{run_id}/demux/done"),
    params:
        infq="analysis/{run_id}/bam2fastq/nomatch_rawlib.basecaller.fastq",
        outfq="analysis/{run_id}/demux/%.fq",
        IDfq="analysis/{run_id}/demux/ID.fq",
        unmatched_fq = "analysis/{run_id}/demux/unmatched.fq",
    message: "Execute the demultiplex step by fastq-multx",
    log:
        stdout="analysis/{run_id}/logs/demux/{run_id}.o",
        stderr="analysis/{run_id}/logs/demux/{run_id}.o",
    shell:
        """
        fastq-multx -m 0 -b -B {input.barcode_info} {params.infq} -o {params.outfq} 2>{log.stderr} 1>{log.stdout}
        rm {params.IDfq}
        rm {params.unmatched_fq} 
        """


rule fastq_scan:
    """
    Run fastq-scan (get fastq stat in CLI)
    """
    input:
        rules.demux.output.touch
    # lambda wildcards: expand("analysis/{run_id}/demux/{sample}.fq",run_id = config['run_id'], sample=get_sample_name()),
    output:
        out_csv="analysis/{run_id}/fastq_scan/{sample}.csv",
        out_summ_csv="analysis/{run_id}/fastq_scan/{sample}.summary.csv",
    message: "Execute fastq-scan for raw data - NOT trimmed data - on {input}",
    log:
        stdout="analysis/{run_id}/logs/fastq-scan/{sample}.o",
        stderr="analysis/{run_id}/logs/fastq-scan/{sample}.e",
    params:
        in_fq="analysis/{run_id}/demux/{sample}.fq",run_id=config['run_id'],sample=get_sample_name(),
    shell:
        """
        cat {params.in_fq} |  fastq-scan  | dasel_linux_386 -r json ".qc_stats" | dasel_linux_386 -r json -w csv > {output.out_csv}
        
        awk -F "," '{{print $2","$6","$3","$9","$13","$10}}' {output.out_csv} | sed 's/\.00*//g' > {output.out_summ_csv}
        """


rule merged_fastq_stat:
    """
    merge fastq-scan outputs
    """
    input:
        # lambda wildcards: expand("analysis/{run_id}/fastq_scan/{sample}.summary.csv",run_id = config['run_id'], sample=get_sample_name()),
        expand("analysis/{run_id}/fastq_scan/{sample}.summary.csv",run_id=config['run_id'],sample=get_sample_name()),
    output:
        touch("analysis/{run_id}/fastq_scan/done"),
    params:
        path="analysis/{run_id}/fastq_scan/",
        summary_file="analysis/{run_id}/fastq_scan/fastq_stat.csv",
    run:
        csv_list = [i for i in glob.glob(params.path + "*summary.csv") if os.path.isfile(i)]
        with open(params.summary_file,'w') as fout:
            fout.write("sample,quality_score_25th,quality_score_median,quality_score_75th,read_length_25th,read_length_median,read_length_75th\n")
            for csv_f in csv_list:
                with open(csv_f) as fin:
                    file_name = os.path.basename(csv_f).split(".")[0]
                    fin.readline()
                    fout.write(f"{file_name},{fin.readline()}")


rule fastqc:
    """
    Run fastqc on raw_data/files.
    """
    input:
        rules.demux.output.touch,
    output:
        html="analysis/{run_id}/fastqc/{sample}_fastqc.html",
        zip="analysis/{run_id}/fastqc/{sample}_fastqc.zip"
    params:
        outdir=os.path.join("analysis/{run_id}/fastqc/"),
        input="analysis/{run_id}/demux/{sample}.fq",run_id=config["run_id"],sample=get_sample_name(),
    message: "Execute fastqc for raw data - NOT trimmed data - on {input}"
    log:
        stdout="analysis/{run_id}/logs/fastqc/{sample}.o",
        stderr="analysis/{run_id}/logs/fastqc/{sample}.e"
    threads: 1
    resources:
        mem_gb=8
    container:
        config["fastqc"]
    # "docker://biocontainers/fastqc:v0.11.9_cv8"
    shell:
        """
        fastqc --outdir {params.outdir} {params.input} 2>{log.stderr} 1>{log.stdout}
        """


rule fastq_screen:
    input:
        rules.demux.output.touch,
    output:
        html="analysis/{run_id}/fastq_screen/{sample}_screen.html",
        txt="analysis/{run_id}/fastq_screen/{sample}_screen.txt",
    log:
        stderr="analysis/{run_id}/logs/fastq_screen/{sample}.e",
        stdout="analysis/{run_id}/logs/fastq_screen/{sample}.o",
    params:
        fastq_screen=config["fastq_screen"]["path"],
        conf=config["fastq_screen"]["conf"],
        bowtie2=config["bowtie2"]["path"],
        input="analysis/{run_id}/demux/{sample}.fq",run_id=config["run_id"],sample=get_sample_name(),
    message: "Execute fastq-screen for checking if samples are contaminated by other random species. - on {input}"
    threads: 2
    container:
        config["fastq_screen"]["sif"]
    # "docker://quay.io/biocontainers/fastq-screen:0.14.0--pl5321hdfd78af_2"
    shell:
        """
        {params.fastq_screen} --aligner bowtie2 \ 
        --bowtie2 {params.bowtie2} \ 
        --outdir analysis/fastq_screen/ \ 
        --threads {threads} \ 
        --conf {params.conf}  {params.input} 2> {log.stderr} 1> {log.stdout}
        """


#deprecated. very few reads left after filtering out low Q reads.
rule fastp:
    """
    Run fastp - this is a quality trimming step.
    """
    input:
        rules.demux.output.touch
    output:
        trimmed_fq="analysis/{run_id}/fastp/{sample}.fq",
        outjson="analysis/{run_id}/fastp/{sample}.fastp.json",
        outhtml="analysis/{run_id}/fastp/{sample}.fastp.html",
    log:
        stdout="logs/{run_id}/fastp/{sample}.o",
        stderr="logs/{run_id}/fastp/{sample}.e",
    threads: 4
    resources:
        mem_gb=8
    params:
        in_fq="analysis/{run_id}/demux/{sample}.fq",run_id=config["run_id"],sample=get_sample_name(),
    container:
        config["fastp"],
    # "docker://quay.io/biocontainers/fastp:0.23.2--h79da9fb_0"
    shell:
        """
        fastp -i {params.in_fq} -o {output.trimmed_fq} -l 200 --json {output.outjson} --html {output.outhtml} --thread 3 --adapter_sequence GTAC
        """

#deprecated
rule cutadapt:
    """
    Run cutadapt to remove short reads/adapter - less than 200bp 
    """
    input:
        rules.demux.output.touch
    output:
        trimmed_fq="analysis/{run_id}/cutadapt/{sample}.fq",
    log:
        stdout="analysis/{run_id}/logs/cutadapt/{sample}.o",
        stderr="analysis/{run_id}/logs/cutadapt/{sample}.e",
    params:
        length="200",
        adapter="GTAC",
        input="analysis/{run_id}/demux/{sample}.fq",run_id=config["run_id"],sample=get_sample_name(),
    message: "Execute cutadapt to trim the adapter sequences as well as filter out reads shorter than 200 bp."
    threads: 4
    resources:
        mem_gb=8
    container:
        config["cutadapt"],
    # "docker://quay.io/biocontainers/cutadapt:4.1--py38hbff2b2d_1"
    shell:
        """
        cutadapt -o {output.trimmed_fq} \
        --minimum-length {params.length} \
        -b {params.adapter} \
        --revcomp {params.input} 1> {log.stdout} 2> {log.stderr}
        """

# no longer needed.
rule pollux:
    """
    Run pollux to correct homopolymer errors.
    """
    input:
        rules.cutadapt.output.trimmed_fq,
    output:
        fqgz="analysis/{run_id}/pollux/{sample}.fq.gz",
        tempfq=temp("analysis/{run_id}/pollux/{sample}.fq.corrected"),
        lowfq=temp("analysis/{run_id}/pollux/{sample}.fq.low"),
    message: "Execute pollux to correct short indel errors on {input}"
    log:
        stderr="analysis/{run_id}/logs/pollux/{sample}.e",
        stdout="analysis/{run_id}/logs/pollux/{sample}.o",
    params:
        outdir="analysis/{run_id}/pollux",
    shell:
        """
        pollux -i {input} -o {params.outdir} -h "true" 2> {log.stderr} 1> {log.stdout}
        
        gzip -c {output.tempfq} > {output.fqgz}
        """


rule bwa:
    """
    Run bwa-mem
    """
    input:
        fq=rules.cutadapt.output.trimmed_fq,
        idx=config['ref_modi']['index'],
    output:
        outsam=temp("analysis/{run_id}/bwamem/{sample}.sam")
    params:
        prefix="{sample}",
    log:
        stdout="logs/{run_id}/bwamem/{sample}.o",
        stderr="logs/{run_id}/bwamem/{sample}.e",
    message: "bwa on {input}"
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

# -M: mark shorter split hits as secondary
# This is optional for Picard compatibility as MarkDuplicates can directly process BWA's alignment, whether or not the alignment marks secondary hits. However, if we want MergeBamAlignment to reassign proper pair alignments, to generate data comparable to that produced by the Broad Genomics Platform, then we must mark secondary alignments [GATK discussion forum].


rule sam2bam:
    input:
        insam=rules.bwa.output.outsam,
    output:
        outbam="analysis/{run_id}/bwamem/{sample}.bam",
        outbai="analysis/{run_id}/bwamem/{sample}.bam.bai",
        idxstat="analysis/{run_id}/bwamem/{sample}.bam.idxstat",
    log:
        stdout="analysis/{run_id}/logs/sam2bam/{sample}.o",
        stderr="analysis/{run_id}/logs/sam2bam/{sample}.e",
    threads: 4
    resources:
        mem_gb=20
    message: "Sam file into bam file - on {input}"
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
        out_stats="analysis/{run_id}/samtools_stats/{sample}.stats",
        out_flagstats="analysis/{run_id}/samtools_stats/{sample}.flagstat.txt",
    threads: 2
    container:
        config["samtools"],
    message: "Caculate alignment stats - on {input}"
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
        output_coverage_tsv="analysis/{run_id}/plot_coverage/{sample}.coverage.tsv",
        tempbam=temp("analysis/{run_id}/samtools{sample}.temp.bam"),
        outbam=("analysis/{run_id}/samtools/{sample}.bam"),
        idxstat="analysis/{run_id}/samtools/{sample}.bam.idxstat",
        bam_idx="analysis/{run_id}/samtools/{sample}.bam.bai",
    log:
        stdout="logs/{run_id}/samtools/{sample}.o",
        stderr="logs/{run_id}/samtools/{sample}.e",
    message: "Filter out supplementary alignment reads and 2nd alignment reads in bam files - on {input.bam}"
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


# Not executed
rule get_align_metrics:
    """
    Run picard to get alignment metrics
    https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering
    """
    input:
        inbam=rules.samtools.output.outbam
    output:
        metrics="analysis/{run_id}/align_metrics/{sample}.align.metrics.txt",
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


# Not executed
# this is an amplicon seq. most of reads have nearly the same genomic coordinate.
rule mark_dups:
    """
    Run picard to mark duplicated reads
    https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering
    """
    input:
        inbam=rules.samtools.output.outbam,
        metrics=rules.get_align_metrics.output.metrics,
    output:
        outbam="analysis/{run_id}/mark_dups/{sample}.bam",
        outbamidx="analysis/{run_id}/mark_dups/{sample}.bam.bai",
    log:
        out="logs/{run_id}/mark_dup/{sample}.o",
        err="logs/{run_id}/mark_dup/{sample}.e"
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
        outvcf="analysis/{run_id}/haplotype_caller/{sample}.vcf",
    params:
        ref_fa=config["ref_modi"]["sequence"],
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
        outvcf="analysis/{run_id}/mutect2/{sample}.vcf",
    params:
        ref_fa=config["ref_modi"]["sequence"],
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
        expand("analysis/{run_id}/fastqc/{sample}_fastqc.html",sample=get_sample_name(),run_id=config["run_id"]),
        expand("analysis/{run_id}/fastqc/{sample}_fastqc.zip",sample=get_sample_name(),run_id=config["run_id"]),
        # expand("analysis/bwamem/{sample}.bam.idxstat", sample=get_sample_name()),
        # expand("analysis/samtools/{sample}.bam.idxstat",sample=get_sample_name()),
        # expand("analysis/fastp/{sample}.fastp.json",sample=get_sample_name()),
        # expand("analysis/samtools/{sample}.flagstats.tsv",sample=get_sample_name()),
        # expand("analysis/align_metrics/{sample}.align.metrics.txt",sample=get_sample_name()),
        # expand("logs/mark_dup/{sample}.o",sample=get_sample_name()),
        expand("analysis/{run_id}/samtools_stats/{sample}.stats",sample=get_sample_name(),run_id=config["run_id"]),
    # expand("analysis/fastq_screen/{sample}_screen.txt",sample=get_sample_name()),
    # expand("analysis/samtools_stats/{sample}.flagstat.txt", sample=get_sample_name()),
    output:
        "analysis/{run_id}/multiqc/multiqc_report.html",
    log:
        stdout="analysis/{run_id}/logs/multiqc/multiqc.o",
        stderr="analysis/{run_id}/logs/multiqc/multiqc.e",
    params:
        outdir="analysis/{run_id}/multiqc/"
    threads: 4
    resources:
        mem_gb=100
    container:
        config["multiqc"],
    # "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0"
    shell:
        """
        multiqc -f {input}  -o {params.outdir} -n multiqc_report.html --cl-config 'max_table_rows: 999999' 2>{log.stdout}
        """


rule bcftools:
    """
    Run bcftools to call variants
    """
    input:
        inbam=rules.samtools.output.outbam
    output:
        outvcf="analysis/{run_id}/bcftools/{sample}.vcf"
    log:
        stdout="analysis/{run_id}/logs/bcftools/{sample}.o",
        stderr="analysis/{run_id}/logs/bcftools/{sample}.e",
    params:
        ref_fa=config['ref_modi']['index'],
    container:
        config['bcftools'],
    # "docker://quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0"
    shell:
        """
        bcftools mpileup -f {params.ref_fa} {input.inbam} | bcftools call -mv -Ov -o {output.outvcf} 2> {log.stdout} 1> {log.stderr}
        """


rule freebayes:
    """
    Run freebayes to call variants
    """
    input:
        inbam=rules.samtools.output.outbam
    output:
        outvcf="analysis/{run_id}/freebayes/{sample}.vcf"

    params:
        ref_fa=config['ref_modi']['index'],
    container:
        config["freebayes"]
    # "docker://quay.io/biocontainers/freebayes:1.3.6--hb089aa1_0"
    shell:
        """
        freebayes -f {params.ref_fa} -C 50 --min-mapping-quality 4 --read-max-mismatch-fraction 1 --read-snp-limit 10 {input.inbam} > {output.outvcf}
        """

#freebayes -f {params.ref_fa} -C 10 --min-mapping-quality 4 --read-max-mismatch-fraction 1 --min-coverage 6 --read-snp-limit 10 {input.inbam} > {output.outvcf}
# this is freebayes options thresholds.
# go to http://194.141.43.170/ion-docs/GUID-D6EA72BF-74FE-4923-B796-E484734B424B.html


rule vcftools:
    """
    Filter vcf files from freebayes
    """
    input:
        rules.freebayes.output.outvcf,
    output:
        out_filt_vcf="analysis/{run_id}/vcftools/{sample}.filt.recode.vcf",
    params:
        minQ=10,
        minGQ=10,
        prefix="analysis/{run_id}/vcftools/{sample}.filt",
    message: "Filter out variants having QUAL less than 10 and GQ less than 10. - on {input}"
    log:
        stdout="analysis/{run_id}/logs/vcftools/{sample}.o",
        stderr="analysis/{run_id}/logs/vcftools/{sample}.e",
    container:
        config["vcftools"],
    shell:
        """
        vcftools --vcf {input} --minQ {params.minQ} --recode --recode-INFO-all --out {params.prefix} --minGQ {params.minGQ} 2> {log.stderr} 1> {log.stdout}
        """

# diff between QUAL and GQ
# see https://gatk.broadinstitute.org/hc/en-us/articles/360035531392-Difference-between-QUAL-and-GQ-annotations-in-germline-variant-calling


rule variant_report:
    """
    VCF file reports - html
    """
    input:
        rules.vcftools.output.out_filt_vcf
    output:
        outhtml="analysis/{run_id}/vcf_html/{sample}.vcf.html",
    message: "Generate human friendly vcf reports on {input}"
    params:
        jar="/usr/local/bin/vcf2table.jar"
    shell:
        """
        java -jar {params.jar} {input} --color --format html > {output.outhtml}
        """


rule deepvariant:
    input:
        inbam=rules.samtools.output.outbam
    output:
        outvcf="analysis/{run_id}/deepvariant/{sample}.vcf",
        outgvcf="analysis/{run_id}/deepvariant/{sample}.gvcf",
        outhtml="analysis/{run_id}/deepvariant/{sample}.visual_report.html"
    message: "Run deepvariant on {input.inbam} - this tool is used for mainly germline variants calling"
    log:
        stdout="analysis/{run_id}/logs/deepvariant/{sample}.o",
        stderr="analysis/{run_id}/logs/deepvariant/{sample}.e",
    params:
        ref_fa=config["ref_modi"]["sequence"],
        deepvariant_sif=config["deepvariant"],
        threads=4,
        region="exon_4_to_10:1-1129",
        ref_fa_path=config["ref_modi"]["path"],
    shell:
        """
        singularity run -B {params.ref_fa_path}:{params.ref_fa_path} {params.deepvariant_sif} /opt/deepvariant/bin/run_deepvariant \
         --model_type=WES \
         --ref="{params.ref_fa}" \
         --reads={input.inbam} \
         --regions "{params.region}" \
         --output_vcf={output.outvcf} \
         --output_gvcf={output.outgvcf} \
         --num_shards={params.threads} 2> {log.stderr} 1> {log.stdout}
        """


rule get_coverage:
    """
    Run David's script
    """
    input:
        inbam=rules.samtools.output.outbam,
    output:
        summaryfile="analysis/{run_id}/get_coverage/{sample}.txt",
        bedfile="analysis/{run_id}/get_coverage/{sample}.tsv",
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
        outpng="analysis/{run_id}/plot_coverage/{sample}.png",
    params:
        prefix="analysis/{run_id}/plot_coverage/{sample}"
    container:
        config["tidyverse"]
    shell:
        """
        Rscript src/coverage.R {input.coverage_tsv} {params.prefix}
        """


#deprecated. after filtering, many reads filtered out.



rule STAR:
    input:
        rules.cutadapt.output.trimmed_fq
    output:
        # see STAR manual for additional output files
        bam="analysis/{run_id}/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai="analysis/{run_id}/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        log_final="analysis/{run_id}/star/{sample}.Log.final.out",
        log="analysis/{run_id}/star/{sample}.Log.out",
        rpg="analysis/{run_id}/star/{sample}.ReadsPerGene.out.tab",
        sj="analysis/{run_id}/star/{sample}.SJ.out.tab",
        g_dir=directory("analysis/{run_id}/star/{sample}._STARgenome"),
        pass1_dir=directory("analysis/{run_id}/star/{sample}._STARpass1"),
    params:
        # path to STAR reference genome index
        index=config["ref_modi"]["star_index"],
        outprefix="analysis/{run_id}/star/{sample}."
    log:
        "logs/{run_id}/star/{sample}.log"
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
        outbam="analysis/{run_id}/star_samtools/{sample}.bam",
        bamidx="analysis/{run_id}/star_samtools/{sample}.bam.bai"
    threads: 4
    log:
        stdout="logs/{run_id}/star_samtools/{sample}.o",
        stderr="logs/{run_id}/star_samtools/{sample}.e",
    container:
        "docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    shell:
        """
        # If you want to keep uniquely-mapped reads + MAPQ > 20 #
        samtools view -h -F 256 -q 20 {input} | samtools sort -O BAM -o {output.outbam} > {log.stdout} 2> {log.stderr}
        
        samtools index -@ {threads} {output.outbam} 
        """
