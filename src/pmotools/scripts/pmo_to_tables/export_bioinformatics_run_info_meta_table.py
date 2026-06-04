#!/usr/bin/env python3
import argparse
import sys


from pmotools.pmo_engine.pmo_exporter import PMOExporter
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def get_parser_export_bioinformatics_run_info_meta_table() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pmotools-python export_bioinformatics_run_info_meta_table",
        description="export the bioinformatics_run_info meta table from a PMO file",
    )
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
    return parser


def parse_args_export_bioinformatics_run_info_meta_table():
    parser = get_parser_export_bioinformatics_run_info_meta_table()
    return parser.parse_args()


def export_bioinformatics_run_info_meta_table():
    args = parse_args_export_bioinformatics_run_info_meta_table()

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
    info_df = PMOExporter.export_bioinformatics_run_info_meta_table(pmo)

    # output
    info_df.to_csv(
        sys.stdout if "STDOUT" == args.output else args.output,
        sep=output_delim,
        index=False,
    )


if __name__ == "__main__":
    export_bioinformatics_run_info_meta_table()
