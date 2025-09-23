#!/usr/bin/env python3
import argparse
import sys


from pmotools.pmo_engine.pmo_processor import PMOProcessor
from pmotools.pmo_engine.pmo_reader import PMOReader
from pmotools.utils.small_utils import Utils


def parse_args_count_specimen_meta():
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
    parser.add_argument(
        "--meta_fields",
        type=str,
        required=True,
        help="the fields to count the subfields of, can supply multiple separated by commas, e.g. --meta_fields collection_country,collection_date",
    )

    return parser.parse_args()


def count_specimen_meta():
    args = parse_args_count_specimen_meta()

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

    # process the meta_fields argument
    meta_fields_toks = args.meta_fields.split(",")

    # read in PMO
    pmo = PMOReader.read_in_pmo(args.file)

    # count sub-fields
    counts_df = PMOProcessor.count_specimen_by_field_value(pmo, meta_fields_toks)

    # write out
    counts_df.to_csv(
        sys.stdout if "STDOUT" == args.output else args.output,
        sep=output_delim,
        index=False,
    )


if __name__ == "__main__":
    count_specimen_meta()
