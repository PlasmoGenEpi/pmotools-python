#!/usr/bin/env python3
import os, argparse, json
import sys
from collections import defaultdict

import pandas as pd

from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_list_experiment_sample_ids_per_specimen_id():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--output', type=str, default="STDOUT", required=False, help='output file')
    parser.add_argument('--delim', default="tab", type=str, required=False, help='the delimiter of the output text file, examples input tab,comma but can also be the actual delimiter')

    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')

    return parser.parse_args()

def list_experiment_sample_ids_per_specimen_id():
    args = parse_args_list_experiment_sample_ids_per_specimen_id()


    # check files
    output_delim, output_extension = Utils.process_delimiter_and_output_extension(args.delim, gzip=args.output.endswith(".gz"))
    args.output = args.output if "STDOUT" == args.output else Utils.appendStrAsNeeded(args.output, output_extension)
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # count fields
    info_df = PMOProcessor.list_experiment_sample_ids_per_specimen_id(pmo)

    # output
    info_df.to_csv(sys.stdout if "STDOUT" == args.output else args.output, sep = output_delim, index=False)

if __name__ == "__main__":
    list_experiment_sample_ids_per_specimen_id()

