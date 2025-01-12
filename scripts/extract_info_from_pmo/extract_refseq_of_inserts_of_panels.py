#!/usr/bin/env python3
import os, argparse, json
import sys
from collections import defaultdict

import pandas as pd

from pmotools.pmo_utils.PMOExtractor import PMOExtractor
from pmotools.pmo_utils.PMOReader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_extract_refseq_of_inserts_of_panels():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--output', type=str, default="STDOUT", required=False, help='output file')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    parser.description = "extract ref_seq of inserts of panels, but if no ref_seq is save in the PMO will just be blank"
    return parser.parse_args()

def extract_refseq_of_inserts_of_panels():
    args = parse_args_extract_refseq_of_inserts_of_panels()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # get panel insert locations
    panel_bed_locs = PMOExtractor.extract_panels_insert_bed_loc(pmo)

    # write
    output_target = sys.stdout if args.output == "STDOUT" else open(args.output, "w")
    with output_target as f:
        f.write("\t".join(["panel_id", "target_id", "ref_seq"]) + "\n")
        for panel_id, bed_locs in panel_bed_locs.items():
            for loc in bed_locs:
                f.write("\t".join([str(panel_id), loc.name, loc.ref_seq]) + "\n")

if __name__ == "__main__":
    extract_refseq_of_inserts_of_panels()

