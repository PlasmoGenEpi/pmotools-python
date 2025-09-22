#!/usr/bin/env python3
import argparse
import sys


from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_list_bioinformatics_run_names():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="PMO file")
    parser.add_argument(
        "--output", type=str, default="STDOUT", required=False, help="output file"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )

    return parser.parse_args()


def list_bioinformatics_run_names():
    args = parse_args_list_bioinformatics_run_names()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # extract all bio run names
    bio_run_names = PMOProcessor.get_bioinformatics_run_names(pmo)

    # write
    output_target = sys.stdout if args.output == "STDOUT" else open(args.output, "w")
    with output_target as f:
        f.write("\n".join(bio_run_names) + "\n")


if __name__ == "__main__":
    list_bioinformatics_run_names()
