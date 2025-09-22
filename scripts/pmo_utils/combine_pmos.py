#!/usr/bin/env python3
import argparse


from pmotools.pmo_engine.pmo_writer import PMOWriter
from pmotools.utils.small_utils import Utils
from pmotools.pmo_engine.pmo_reader import PMOReader


def parse_args_combine_pmos():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pmo_files",
        type=str,
        required=True,
        help="a list of PMO files to combine into 1 PMO file, must be from same amplicon panel",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Output new combined PMO file"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )

    return parser.parse_args()


def combine_pmos():
    args = parse_args_combine_pmos()

    # set up output
    args.output = PMOWriter.add_pmo_extension_as_needed(
        args.output, args.output.endswith(".gz")
    )
    Utils.outputfile_check(args.output, args.overwrite)

    # check if at least 2 PMO files supplied
    pmo_files_list = Utils.parse_delimited_input_or_file(args.pmo_files, ",")
    if len(pmo_files_list) < 2:
        raise Exception(
            "Only supplied "
            + str(len(pmo_files_list))
            + " but multiple PMO files were expected"
        )

    # read in the PMOs
    pmos = PMOReader.read_in_pmos(pmo_files_list)

    # combine PMOs
    pmo_out = PMOReader.combine_multiple_pmos(pmos)

    # write
    PMOWriter.write_out_pmo(pmo_out, args.output, args.overwrite)


if __name__ == "__main__":
    combine_pmos()
