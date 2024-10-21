#!/usr/bin/env python3
import os, argparse, json
import sys
from collections import defaultdict

import pandas as pd

from pmotools.extract_from_pmo.PMOExtractor import PMOExtractor
from pmotools.extract_from_pmo.PMOReader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_count_targets_per_sample():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--output', type=str, default="STDOUT", required=False, help='output file')
    parser.add_argument('--delim', default="tab", type=str, required=False, help='the delimiter of the output text file, examples input tab,comma but can also be the actual delimiter')
    parser.add_argument('--overwrite', action='store_true', help='If output file exists, overwrite it')
    parser.add_argument('--read_count_minimum', default=0.0, type=float, required=False, help='the minimum read count (inclusive) to be counted as covered by sample')

    return parser.parse_args()


def count_targets_per_sample():
    args = parse_args_count_targets_per_sample()

    output_delim, output_extension = Utils.process_delimiter_and_output_extension(args.delim)

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # count
    counts_df = PMOExtractor.count_targets_per_sample(pmo, args.read_count_minimum)

    if "STDOUT" == args.output:
        counts_df.to_csv(sys.stdout, sep=output_delim, index=False)
    else:
        counts_df.to_csv(args.output, sep=output_delim, index=False)


if __name__ == "__main__":
    count_targets_per_sample()

