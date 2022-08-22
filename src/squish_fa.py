#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/9/22

@author: jaykim
"""
import argparse
import os.path
from dataclasses import dataclass
from typing import Dict, List

from extract_exons import exists


@dataclass
class Args:
    ref_fa: str
    out_fa: str


def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description="combine multiple ids/seqs from a fasta file into one id/seq",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i",
                        "--ref_fa",
                        help="fasta file",
                        metavar="file",
                        type=str,
                        required=True
                        )

    parser.add_argument("-o",
                        "--out_fa",
                        help="fasta file",
                        metavar="file",
                        type=str,
                        required=True
                        )
    args = parser.parse_args()

    return Args(args.ref_fa, args.out_fa)


def read_fa(fa: str) -> Dict:
    result: Dict[str, str] = dict()

    exists(fa)

    with open(fa, "rt") as fin:
        header: List[str] = list()
        seq: List[str] = list()
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                head: str = line.replace(">", "") + "|"
                header.append(head)
            else:
                seq.append(line)

    result["".join(header)] = "".join(seq)
    # print(result)
    return result


def write_fa(fa_dict: Dict[str, str], fa: str, out_fa_prefix: str) -> None:
    exists(fa)
    output_dir: str = os.path.dirname(fa)
    tmp_fa: str = os.path.basename(fa)
    split_fa_name = os.path.splitext(tmp_fa)
    # fa_prefix: str = split_fa_name[0]
    fa_exten: str = split_fa_name[1]

    output_fa: str = os.path.join(output_dir, out_fa_prefix + fa_exten)

    with open(output_fa, "wt") as fout:
        for key in fa_dict:
            # fout.write(">" + key + "\n")
            fout.write(">" + "exon_4_5_6_7" + "\n")
            fout.write(fa_dict.get(key) + "\n")


def main():
    args = get_args()
    ref_fa: str = args.ref_fa
    out_fa_prefix: str = args.out_fa
    fa_dict: Dict[str, str] = read_fa(ref_fa)
    write_fa(fa_dict, ref_fa, out_fa_prefix)


if __name__ == '__main__':
    main()
