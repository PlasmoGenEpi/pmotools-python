#!/usr/bin/env python3
import argparse
import json
import pandas as pd
from pmotools.pmo_builder.mhap_table_to_pmo import mhap_table_to_pmo
from pmotools.utils.small_utils import Utils


def parse_args_text_meta_to_json_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Input text file path")
    parser.add_argument(
        "--delim",
        default="tab",
        type=str,
        required=False,
        help="the delimiter of the input text file, examples tab,comma",
    )
    parser.add_argument(
        "--index_col_name",
        type=str,
        required=False,
        help="by default output is a list, if an index column name is supplied it will be a dict with this column as index",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Output json file path"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If output file exists, overwrite it"
    )
    return parser.parse_args()


def text_meta_to_json_meta():
    args = parse_args_text_meta_to_json_meta()
    if args.delim == "tab":
        args.delim = "\t"
    elif args.delim == "comma":
        args.delim = ","
    args.output = Utils.appendStrAsNeeded(args.output, ".json")
    # check if input file exists and if output file exists check if --overwrite flag is set
    Utils.inputOutputFileCheckFromArgParse(args)
    index_col_name = None
    if args.index_col_name is not None:
        index_col_name = args.index_col_name
    contents = pd.read_csv(args.file, sep=args.delim, index_col=index_col_name)

    contents_json = mhap_table_to_pmo(contents, args.index_col_name)

    json.dump(contents_json, open(args.output, "w"), indent=4)


if __name__ == "__main__":
    text_meta_to_json_meta()
