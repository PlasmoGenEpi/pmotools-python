#!/usr/bin/env python3
import argparse
import json
import pandas as pd
from pmotools.pmo_builder.metatable_to_pmo import pandas_table_to_json
from pmotools.utils.small_utils import Utils


def parse_args_excel_meta_to_json_meta():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True, help="Input excel file path")
    parser.add_argument(
        "--sheet",
        type=str,
        required=False,
        help="The sheet to convert, if none provided will default to first sheet",
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


def excel_meta_to_json_meta():
    args = parse_args_excel_meta_to_json_meta()
    args.output = Utils.appendStrAsNeeded(args.output, ".json")
    # check if input file exists and if output file exists check if --overwrite flag is set
    Utils.inputOutputFileCheckFromArgParse(args)

    sheet = 1
    index_col_name = None
    if args.sheet is not None:
        sheet = args.sheet
    if args.index_col_name is not None:
        index_col_name = args.index_col_name
    contents = pd.read_excel(args.file, sheet_name=sheet, index_col=index_col_name)

    contents_json = pandas_table_to_json(contents, args.index_col_name)

    json.dump(contents_json, open(args.output, "w"), indent=4)


if __name__ == "__main__":
    excel_meta_to_json_meta()
