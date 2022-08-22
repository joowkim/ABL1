#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/8/22

@author: jaykim
"""
import argparse
import json
import os.path
import shutil
import subprocess
from dataclasses import dataclass
from typing import List

import pandas as pd


@dataclass
class Args:
    trans_id: str
    ref_gene_gtf: str
    ref_fa_file: str
    exon_nums: List[str]


def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description="Get exon info based on transcript ID",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-g",
                        "--gtf",
                        help="GTF file",
                        metavar="file",
                        type=str,
                        default="/media/hd1/Repository/CFTR/usrdef/RefGene/refGene_05_14_2020.txt"
                        )

    parser.add_argument("-t",
                        "--trans_id",
                        help="Transcript ID",
                        type=str,
                        )

    parser.add_argument("-rf",
                        "--ref_fa",
                        help="Reference fasta file",
                        type=str,
                        metavar="file",
                        default="/media/hd1/Repository/Cancer/dbs/genome/assembly/ucsc/hg19/hg19.fa"
                        )

    parser.add_argument("--exon_nums",
                        help="Exon numbers you want to extract",
                        type=str,
                        metavar="exon_nums",
                        nargs="+"
                        )

    args = parser.parse_args()

    return Args(args.trans_id, args.gtf, args.ref_fa, args.exon_nums)


def exists(string: str) -> None:
    assert os.path.exists(string), f"{string} is not found!"


def get_gtf_colname(ref_gtf_col_json: str) -> List:
    """

    :param ref_gtf_col_json:
    :return: "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"
    """
    exists(ref_gtf_col_json)
    data = None
    with open(ref_gtf_col_json) as fin:
        data = json.load(fin)
    return list(data.values())[0]


def get_exon_info(trans_id: str, ref_gene_gtf: str) -> pd.DataFrame:
    """

    :param trans_id:
    :return:
     """
    ## this is fur testing
    # ref_gene_gtf: str = "/media/hd1/Repository/CFTR/usrdef/RefGene/refGene_05_14_2020.txt"
    # trans_id = "NM_005157"

    exists(ref_gene_gtf)
    # assert os.path.isfile(ref_gene_gtf), f"{ref_gene_gtf} is not found!"
    tmp_list: List = list()

    with open(ref_gene_gtf, "rt") as fin:
        for r in fin:
            name, transcript, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, \
            exonCount, exonStarts, exonEnds, exonFrames, gene, str1, str2, str3 = r.split("\t")
            # if not transcript.split(".")[0] == trans_id:
            #     continue
            if transcript.split(".")[0] == trans_id:
                exon_starts = list(map(int, exonStarts.strip(",").split(",")))
                exon_ends = list(map(int, exonEnds.strip(",").split(",")))
                # conversion from 0-based to 1-based
                for idx in range(len(exon_starts)):
                    no = idx + 1
                    # exon_starts[idx] += 1
                    exonstart = exon_starts[idx]
                    # print(exonstart)
                    exonend = exon_ends[idx]
                    # exon_regions = zip(exon_starts, exon_ends)
                    tmp_list.append([gene, chrom, strand, transcript, exonstart, exonend, str(no)])

    df = pd.DataFrame(tmp_list, columns=["Gene_ID", "Chrm", "Strand", "Transcript", "ExonStart", "ExonEnd",
                                         "Exon_Number"]).set_index("Gene_ID")
    return df


def select_exons(exon_info_df: pd.DataFrame, exon_numbers: List[str]) -> pd.DataFrame:
    return exon_info_df.loc[exon_info_df["Exon_Number"].isin(exon_numbers)]


def make_bed_file(exon_info_df: pd.DataFrame, trans_id: str, out_bed: str) -> None:
    """
    use bedtools get sequences from fa file. a bed file is required.
    :return:
    """
    ## example bed file
    # chr1 5 10

    result_df: pd.DataFrame = exon_info_df[["Chrm", "ExonStart", "ExonEnd", "Exon_Number"]]
    # out_bed_file_path: str = os.path.join("output", out_bed)
    result_df.to_csv(out_bed, sep="\t", index=False, header=False)
    print(f"generate {out_bed}!")


def get_sequences(bed_file: str, ref_fa: str, out_fa: str) -> None:
    """
    bedtools getfasta -fi test.fa -bed test.bed -fo test.fa.out
    :param bed_file:
    :param ref_fa:
    :param out_fa:
    :return:
    """
    exists(bed_file)
    exists(ref_fa)

    bedtools: str = "bedtools"
    sub_out = subprocess.run(f"{bedtools}")

    ## check if bedtools is not installed.
    assert sub_out.returncode == 0, "bedtools is not found!"

    # os.makedirs("out", exist_ok=True)

    run_out = subprocess.run(["bedtools", "getfasta", "-fi", f"{ref_fa}", "-bed", f"{bed_file}", "-fo", f"{out_fa}"])

    if run_out.returncode == 0:
        print("generating custom_fa is done!")
    else:
        exists(out_fa)
        print("generating custom_fa is NOT done!")
        print("this is the commend executed")
        print(run_out.args)
        print()
        print("stout")
        print(run_out.stdout)
        print()
        print("stderr")
        print(run_out.stderr)


def amend_out_fa(out_fa: str) -> None:
    """

    :param out_fa:
    :return:
    """
    exists(out_fa)

    bak_file: str = f"{out_fa}.bak"
    shutil.copy(out_fa, bak_file)

    with open(bak_file, "rt") as fin, open(out_fa, "wt") as fout:
        for line in fin:
            if line.startswith(">"):
                # line looks like this
                # >chr9:133738149-133738422
                chrm: str = line.split(":")[0]
                tmp_list: List[str] = line.split(":")[1].split("-")
                pos_start: int = int(tmp_list[0])
                pos_end: int = int(tmp_list[1])

                header: str = f"{chrm}:{pos_start + 1}-{pos_end}\n"
                fout.write(header)
            else:
                fout.write(line)


def main():
    # trans_id = "NM_005157"
    # ref_gene_gtf = "/media/hd1/Repository/CFTR/usrdef/RefGene/refGene_05_14_2020.txt"
    # exon_nums = ["4", "5", "6", "7"]
    args = get_args()
    trans_id: str = args.trans_id
    gtf_file: str = args.ref_gene_gtf
    ref_fa: str = args.ref_fa_file
    exon_nums: List[str] = args.exon_nums
    exon_info_df: pd.DataFrame = get_exon_info(trans_id=trans_id, ref_gene_gtf=gtf_file)
    print(exon_info_df)
    select_exons_df: pd.DataFrame = select_exons(exon_info_df=exon_info_df, exon_numbers=exon_nums)
    print(select_exons_df)

    os.makedirs("output", exist_ok=True)
    out_bed_path = os.path.join("output", f"{trans_id}.bed")
    out_fa: str = os.path.join("output", f"{trans_id}.fa")
    make_bed_file(exon_info_df=select_exons_df, trans_id=trans_id, out_bed=out_bed_path)
    get_sequences(out_bed_path, ref_fa, out_fa)
    amend_out_fa(out_fa)


if __name__ == "__main__":
    # trans_id = "NM_005157"
    main()
