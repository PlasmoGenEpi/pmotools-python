#!/usr/bin/env python3
import argparse
import sys


from pmotools.pmo_engine.pmo_exporter import PMOExporter
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_export_project_info_meta_table():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="PMO file")
    parser.add_argument(
        "--output", type=str, default="STDOUT", required=False, help="output file"
    )
    parser.add_argument(
        "--delim",
        default="tab",
        type=str,
        required=False,
        help="the delimiter of the output text file, examples input tab,comma but can also be the actual delimiter",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )

    return parser.parse_args()


def export_project_info_meta_table():
    args = parse_args_export_project_info_meta_table()

    # check files
    output_delim, output_extension = Utils.process_delimiter_and_output_extension(
        args.delim, gzip=args.output.endswith(".gz")
    )
    args.output = (
        args.output
        if "STDOUT" == args.output
        else Utils.appendStrAsNeeded(args.output, output_extension)
    )
    Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # count fields
    info_df = PMOExporter.export_project_info_meta_table(pmo)

    # output
    info_df.to_csv(
        sys.stdout if "STDOUT" == args.output else args.output,
        sep=output_delim,
        index=False,
    )


if __name__ == "__main__":
    export_project_info_meta_table()
