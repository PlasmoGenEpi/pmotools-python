#!/usr/bin/env python3
import argparse
import sys

import pandas as pd

from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_export_population_names():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meta_fields", type=str, required=True, help="meta fields")
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
    parser.add_argument(
        "--export_by_specimen_id",
        action="store_true",
        help="export by specimen name instead of the default is export by library_sample_name",
    )

    return parser.parse_args()


def export_population_names():
    args = parse_args_export_population_names()

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
    if args.export_by_specimen_id:
        # get population
        population_map = PMOProcessor.export_population_info_for_specimen_name(
            pmo, args.meta_fields.split(",")
        )

        population_map_df = pd.DataFrame(
            list(population_map.items()), columns=["library_sample_name", "population"]
        )
        # output
        population_map_df.to_csv(
            sys.stdout if "STDOUT" == args.output else args.output,
            sep=output_delim,
            index=False,
        )
    else:
        # get population
        population_map = PMOProcessor.export_population_info_for_library_sample_name(
            pmo, args.meta_fields.split(",")
        )

        population_map_df = pd.DataFrame(
            list(population_map.items()), columns=["library_sample_name", "population"]
        )
        # output
        population_map_df.to_csv(
            sys.stdout if "STDOUT" == args.output else args.output,
            sep=output_delim,
            index=False,
        )


if __name__ == "__main__":
    export_population_names()
