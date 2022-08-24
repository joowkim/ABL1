#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/12/22

@author: David
"""
import argparse
import os
import subprocess
from dataclasses import dataclass

import numpy
import pandas


@dataclass
class Args:
    bed: str
    bam: str
    min_read_over_lap: float
    min_amplicon_cov: int
    min_sample_cov: int
    min_read_uniform: float
    output: str
    summ_file: str


def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description="Get exon info based on transcript ID",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",
                        "--bed",
                        help="amplicons_target_region_bed_file",
                        metavar="bed file",
                        type=str,
                        default="configs/exon_4_5_6_7.bed",
                        dest="bed",
                        )

    parser.add_argument("-b",
                        "--bam",
                        help="bam file",
                        type=str,
                        required=True,
                        dest="bam"
                        )

    parser.add_argument("-min_read_over_lap",
                        help="min_read_over_lap",
                        type=float,
                        default=0.5,
                        dest="min_read_over_lap"
                        )

    parser.add_argument("-min_amplicon_cov",
                        help="min amplicon coverage",
                        type=int,
                        default=300,
                        dest="min_amplicon_cov",
                        )

    parser.add_argument("-min_sample_cov",
                        help="min sample coverage",
                        type=int,
                        default=100,
                        dest="min_sample_cov",
                        )

    parser.add_argument("-min_read_uniform",
                        help="min read uniform",
                        type=float,
                        default=0.2,
                        dest="min_read_uniform",
                        )

    parser.add_argument("-o"
                        "--output_coverage_file",
                        help="output_coverage_file",
                        type=str,
                        dest="output",
                        required=True
                        )

    parser.add_argument("-s"
                        "--summary_file",
                        help="summary file",
                        type=str,
                        dest="summ_file",
                        required=True,
                        )

    args = parser.parse_args()

    return Args(args.bed,
                args.bam,
                args.min_read_over_lap,
                args.min_amplicon_cov,
                args.min_sample_cov,
                args.min_read_uniform,
                args.output,
                args.summ_file,
                )


def onTargetReads(bam_file, amplicons_target_region_bed_file):
    cmd = '''samtools view -F 0x40 -c %s''' % (bam_file)
    mapped_reads = 'na'
    try:
        mapped_reads = int([x for x in subprocess.getoutput(cmd).split('\n') if x.strip() != ''][0])
    except:
        pass
    cmd = '''samtools view -L %s -c %s''' % (amplicons_target_region_bed_file, bam_file)
    mapped_to_target = int([x for x in subprocess.getoutput(cmd).split('\n') if x.strip() != ''][0])
    on_target = 'na'
    try:
        on_target = float(mapped_to_target) * 100 / mapped_reads
    except:
        pass
    return on_target, mapped_reads


def get_coverage(amplicons_target_region_bed_file,
                 bam_file,
                 minReadoverlap,
                 minAmpliconCov,
                 minSampleCov,
                 minReadUniform,
                 output_coverage_file_tsv,
                 summ_file):
    assert os.path.isfile(bam_file), f"{bam_file} is not found!"
    assert os.path.isfile(amplicons_target_region_bed_file), f"{amplicons_target_region_bed_file} is not found!"

    cmd = f'''bedtools coverage -a {amplicons_target_region_bed_file} -b {bam_file} -counts -f {minReadoverlap} > {output_coverage_file_tsv}'''

    subprocess.call(cmd, shell=True)
    # Strg = '\\t'.join(
    #     ['Temp' + str(ix) if not ix in [3, 8] else ('target' if ix == 3 else 'coverage' + '\\n') for ix in
    #      range(9)])

    Strg = ["chrm", "pos_start", "pos_end", "target", "coverage"]

    # cmd = """sed -i '1s/^/%s/' %s""" % (Strg, output_coverage_file_tsv)
    # cmd = f"cat echo {Strg} {output_coverage_file_tsv}"
    # subprocess.call(cmd, shell=True)
    dfcov = pandas.read_csv(output_coverage_file_tsv, sep='\t', names=Strg, header=None)

    mean_depth = numpy.mean(list(dfcov['coverage']))
    QC = 'Pass' if float(mean_depth) >= minSampleCov else 'Fail'

    dfcov['amplicon:coverage'] = dfcov['target'].astype(str) + ':' + dfcov['coverage'].astype(str)
    dfcov['qc'] = dfcov['amplicon:coverage'].apply(
        lambda x: 'Pass' if int(x.split(':')[1]) >= minAmpliconCov else 'Fail')
    dfcov = dfcov.drop(columns=['amplicon:coverage'])
    df = dfcov.loc[dfcov['coverage'] >= minReadUniform * mean_depth]
    Uniformity = len(df) * 100.0 / len(dfcov)

    OnTarget, mapped_reads = onTargetReads(bam_file, amplicons_target_region_bed_file)

    Amplicons = dfcov.drop(columns=[x for x in dfcov.columns if x.startswith('Temp')])
    Columns = ['mapped_reads', 'on_target', 'mean_depth', 'uniformity', 'qc']
    SummaryCov = [int(mapped_reads), round(OnTarget, 1), int(mean_depth), round(Uniformity, 1), QC]
    SummaryCov = pandas.DataFrame([SummaryCov], columns=Columns)
    SummaryCov.to_csv(summ_file, index=False, sep="\t")
    print(SummaryCov)


if __name__ == "__main__":
    # minAmpliconCov = 300
    # minSampleCov = 100
    # minReadUniform = 0.2
    # # minReadoverlap = 0.6 -> the forward reads start in the middle of exon_4... So if I keep using minReadoverlap 0.6, I get 0 reads/coverage for exon 4
    # minReadoverlap = 0.5
    # amplicons_target_region_bed_file = 'configs/exon_4_5_6_7.bed'
    # bam_file = 'before_gatk_analysis/03.samtools/ABL1.bam'
    # output_coverage_file_tsv = 'test_data/out.bed'

    args = get_args()
    minAmpliconCov: int = args.min_amplicon_cov
    minSampleCov: int = args.min_sample_cov
    minReadUniform: float = args.min_read_uniform
    # minReadoverlap = 0.6 -> the forward reads start in the middle of exon_4... So if I keep using minReadoverlap 0.6, I get 0 reads/coverage for exon 4
    minReadoverlap: float = args.min_read_over_lap
    amplicons_target_region_bed_file: str = args.bed
    bam_file: str = args.bam
    output_coverage_file_tsv:str = args.output
    summ_file:str = args.summ_file

    get_coverage(amplicons_target_region_bed_file,
                 bam_file,
                 minReadoverlap,
                 minAmpliconCov,
                 minSampleCov,
                 minReadUniform,
                 output_coverage_file_tsv,
                 summ_file,
                 )
