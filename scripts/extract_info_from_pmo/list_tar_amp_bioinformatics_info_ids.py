#!/usr/bin/env python3
import os, argparse, json
import sys
from collections import defaultdict

import pandas as pd

from pmotools.pmo_utils.PMOExtractor import PMOExtractor
from pmotools.pmo_utils.PMOReader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_list_tar_amp_bioinformatics_info_ids():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required=True, help='PMO file')
    parser.add_argument('--output', type=str, default="STDOUT", required=False, help='output file')
    parser.add_argument('--overwrite', action = 'store_true', help='If output file exists, overwrite it')

    return parser.parse_args()

def list_tar_amp_bioinformatics_info_ids():
    args = parse_args_list_tar_amp_bioinformatics_info_ids()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # extract all taramp_bioinformatics_ids
    bioids = pmo["taramp_bioinformatics_infos"].keys()

    # write
    if "STDOUT" == args.output:
        sys.stdout.write("\n".join(bioids) + "\n")
    else:
        with open(args.output, "w") as f:
            f.write("\n".join(bioids) + "\n")




if __name__ == "__main__":
    list_tar_amp_bioinformatics_info_ids()

