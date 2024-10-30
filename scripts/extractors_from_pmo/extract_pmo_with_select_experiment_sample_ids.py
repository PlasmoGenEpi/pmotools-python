#!/usr/bin/env python3
import os, argparse, json
import sys
from collections import defaultdict

import pandas as pd

from pmotools.pmo_utils.PMOExtractor import PMOExtractor
from pmotools.pmo_utils.PMOReader import PMOReader
from pmotools.pmo_utils.PMOWriter import PMOWriter
from pmotools.utils.small_utils import Utils


def parse_args_extract_pmo_with_select_experiment_sample_ids():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--output', type=str, required=True, help='Output json file path')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')
    parser.add_argument('--verbose', action = 'store_true', help='write out various messages about extraction')
    parser.add_argument('--experiment_sample_ids', type=str, required=True, help='Can either comma separated experiment_sample_ids, or a plain text file where each line is a experiment_sample_id')
    return parser.parse_args()

def extract_pmo_with_select_experiment_sample_ids():
    args = parse_args_extract_pmo_with_select_experiment_sample_ids()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # parse specimen ids
    all_experiment_sample_ids = Utils.parse_delimited_input_or_file(args.experiment_sample_ids)

    # read in pmo
    pmo = PMOReader.read_in_pmo(args.file)

    # extract
    pmo_out = PMOExtractor.extract_from_pmo_select_experiment_sample_ids(pmo, all_experiment_sample_ids)

    # write out the extracted
    args.output = PMOWriter.add_pmo_extension_as_needed(args.output, args.file.endswith('.gz') or args.output.endswith(".gz"))
    PMOWriter.write_out_pmo(pmo_out, args.output, args.overwrite)



if __name__ == "__main__":
    extract_pmo_with_select_experiment_sample_ids()

