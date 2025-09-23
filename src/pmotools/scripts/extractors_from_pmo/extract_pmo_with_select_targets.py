#!/usr/bin/env python3
import argparse


from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.pmo_engine.pmo_writer import PMOWriter
from pmotools.utils.small_utils import Utils


def parse_args_extract_pmo_with_select_targets():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="PMO file")
    parser.add_argument(
        "--output", type=str, required=True, help="Output json file path"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="write out various messages about extraction",
    )
    parser.add_argument(
        "--targets",
        type=str,
        required=True,
        help="Can either comma separated target_namess, or a plain text file where each line is a target_namess",
    )
    return parser.parse_args()


def extract_pmo_with_select_targets():
    args = parse_args_extract_pmo_with_select_targets()

    # check files
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # parse target names
    all_target_names = Utils.parse_delimited_input_or_file(args.targets)

    # read in pmo
    pmo = PMOReader.read_in_pmo(args.file)

    # extract
    pmo_out = PMOProcessor.filter_pmo_by_target_names(pmo, all_target_names)

    # write out the extracted
    args.output = PMOWriter.add_pmo_extension_as_needed(
        args.output, args.file.endswith(".gz") or args.output.endswith(".gz")
    )
    PMOWriter.write_out_pmo(pmo_out, args.output, args.overwrite)


if __name__ == "__main__":
    extract_pmo_with_select_targets()
