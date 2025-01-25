#!/usr/bin/env python3
import os, argparse, json
import sys
from collections import defaultdict

import pandas as pd

from pmotools.pmo_utils.PMOExtractor import PMOExtractor
from pmotools.pmo_utils.PMOReader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_extract_insert_of_panels():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--output', type=str, default="STDOUT", required=False, help='output file')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    parser.add_argument('--add_ref_seqs', action = 'store_true', help='add ref seqs to the output as ref_seq')

    return parser.parse_args()

def extract_insert_of_panels():
    args = parse_args_extract_insert_of_panels()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # get panel insert locations
    panel_bed_locs = PMOExtractor.extract_panels_insert_bed_loc(pmo)

    # write
    output_target = sys.stdout if args.output == "STDOUT" else open(args.output, "w")
    with output_target as f:
        f.write("\t".join(["#chrom", "start", "end", "target_id", "length", "strand", "extra_info"]))
        if args.add_ref_seqs:
            f.write("\tref_seq")
        f.write("\n")
        for panel_id, bed_locs in panel_bed_locs.items():
            for loc in bed_locs:
                f.write("\t".join([loc.chrom, str(loc.start), str(loc.end), loc.name, str(loc.score), loc.strand, loc.extra_info]))
                if args.add_ref_seqs:
                    f.write("\t" + str(loc.ref_seq))
                f.write("\n")

if __name__ == "__main__":
    extract_insert_of_panels()
